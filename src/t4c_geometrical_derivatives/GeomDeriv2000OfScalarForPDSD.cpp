#include "GeomDeriv2000OfScalarForPDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_pdsd_0(CSimdArray<double>& buffer_2000_pdsd,
                     const CSimdArray<double>& buffer_pdsd,
                     const CSimdArray<double>& buffer_fdsd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_pdsd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdsd

    auto g_x_xx_0_xx = buffer_pdsd[0];

    auto g_x_xx_0_xy = buffer_pdsd[1];

    auto g_x_xx_0_xz = buffer_pdsd[2];

    auto g_x_xx_0_yy = buffer_pdsd[3];

    auto g_x_xx_0_yz = buffer_pdsd[4];

    auto g_x_xx_0_zz = buffer_pdsd[5];

    auto g_x_xy_0_xx = buffer_pdsd[6];

    auto g_x_xy_0_xy = buffer_pdsd[7];

    auto g_x_xy_0_xz = buffer_pdsd[8];

    auto g_x_xy_0_yy = buffer_pdsd[9];

    auto g_x_xy_0_yz = buffer_pdsd[10];

    auto g_x_xy_0_zz = buffer_pdsd[11];

    auto g_x_xz_0_xx = buffer_pdsd[12];

    auto g_x_xz_0_xy = buffer_pdsd[13];

    auto g_x_xz_0_xz = buffer_pdsd[14];

    auto g_x_xz_0_yy = buffer_pdsd[15];

    auto g_x_xz_0_yz = buffer_pdsd[16];

    auto g_x_xz_0_zz = buffer_pdsd[17];

    auto g_x_yy_0_xx = buffer_pdsd[18];

    auto g_x_yy_0_xy = buffer_pdsd[19];

    auto g_x_yy_0_xz = buffer_pdsd[20];

    auto g_x_yy_0_yy = buffer_pdsd[21];

    auto g_x_yy_0_yz = buffer_pdsd[22];

    auto g_x_yy_0_zz = buffer_pdsd[23];

    auto g_x_yz_0_xx = buffer_pdsd[24];

    auto g_x_yz_0_xy = buffer_pdsd[25];

    auto g_x_yz_0_xz = buffer_pdsd[26];

    auto g_x_yz_0_yy = buffer_pdsd[27];

    auto g_x_yz_0_yz = buffer_pdsd[28];

    auto g_x_yz_0_zz = buffer_pdsd[29];

    auto g_x_zz_0_xx = buffer_pdsd[30];

    auto g_x_zz_0_xy = buffer_pdsd[31];

    auto g_x_zz_0_xz = buffer_pdsd[32];

    auto g_x_zz_0_yy = buffer_pdsd[33];

    auto g_x_zz_0_yz = buffer_pdsd[34];

    auto g_x_zz_0_zz = buffer_pdsd[35];

    auto g_y_xx_0_xx = buffer_pdsd[36];

    auto g_y_xx_0_xy = buffer_pdsd[37];

    auto g_y_xx_0_xz = buffer_pdsd[38];

    auto g_y_xx_0_yy = buffer_pdsd[39];

    auto g_y_xx_0_yz = buffer_pdsd[40];

    auto g_y_xx_0_zz = buffer_pdsd[41];

    auto g_y_xy_0_xx = buffer_pdsd[42];

    auto g_y_xy_0_xy = buffer_pdsd[43];

    auto g_y_xy_0_xz = buffer_pdsd[44];

    auto g_y_xy_0_yy = buffer_pdsd[45];

    auto g_y_xy_0_yz = buffer_pdsd[46];

    auto g_y_xy_0_zz = buffer_pdsd[47];

    auto g_y_xz_0_xx = buffer_pdsd[48];

    auto g_y_xz_0_xy = buffer_pdsd[49];

    auto g_y_xz_0_xz = buffer_pdsd[50];

    auto g_y_xz_0_yy = buffer_pdsd[51];

    auto g_y_xz_0_yz = buffer_pdsd[52];

    auto g_y_xz_0_zz = buffer_pdsd[53];

    auto g_y_yy_0_xx = buffer_pdsd[54];

    auto g_y_yy_0_xy = buffer_pdsd[55];

    auto g_y_yy_0_xz = buffer_pdsd[56];

    auto g_y_yy_0_yy = buffer_pdsd[57];

    auto g_y_yy_0_yz = buffer_pdsd[58];

    auto g_y_yy_0_zz = buffer_pdsd[59];

    auto g_y_yz_0_xx = buffer_pdsd[60];

    auto g_y_yz_0_xy = buffer_pdsd[61];

    auto g_y_yz_0_xz = buffer_pdsd[62];

    auto g_y_yz_0_yy = buffer_pdsd[63];

    auto g_y_yz_0_yz = buffer_pdsd[64];

    auto g_y_yz_0_zz = buffer_pdsd[65];

    auto g_y_zz_0_xx = buffer_pdsd[66];

    auto g_y_zz_0_xy = buffer_pdsd[67];

    auto g_y_zz_0_xz = buffer_pdsd[68];

    auto g_y_zz_0_yy = buffer_pdsd[69];

    auto g_y_zz_0_yz = buffer_pdsd[70];

    auto g_y_zz_0_zz = buffer_pdsd[71];

    auto g_z_xx_0_xx = buffer_pdsd[72];

    auto g_z_xx_0_xy = buffer_pdsd[73];

    auto g_z_xx_0_xz = buffer_pdsd[74];

    auto g_z_xx_0_yy = buffer_pdsd[75];

    auto g_z_xx_0_yz = buffer_pdsd[76];

    auto g_z_xx_0_zz = buffer_pdsd[77];

    auto g_z_xy_0_xx = buffer_pdsd[78];

    auto g_z_xy_0_xy = buffer_pdsd[79];

    auto g_z_xy_0_xz = buffer_pdsd[80];

    auto g_z_xy_0_yy = buffer_pdsd[81];

    auto g_z_xy_0_yz = buffer_pdsd[82];

    auto g_z_xy_0_zz = buffer_pdsd[83];

    auto g_z_xz_0_xx = buffer_pdsd[84];

    auto g_z_xz_0_xy = buffer_pdsd[85];

    auto g_z_xz_0_xz = buffer_pdsd[86];

    auto g_z_xz_0_yy = buffer_pdsd[87];

    auto g_z_xz_0_yz = buffer_pdsd[88];

    auto g_z_xz_0_zz = buffer_pdsd[89];

    auto g_z_yy_0_xx = buffer_pdsd[90];

    auto g_z_yy_0_xy = buffer_pdsd[91];

    auto g_z_yy_0_xz = buffer_pdsd[92];

    auto g_z_yy_0_yy = buffer_pdsd[93];

    auto g_z_yy_0_yz = buffer_pdsd[94];

    auto g_z_yy_0_zz = buffer_pdsd[95];

    auto g_z_yz_0_xx = buffer_pdsd[96];

    auto g_z_yz_0_xy = buffer_pdsd[97];

    auto g_z_yz_0_xz = buffer_pdsd[98];

    auto g_z_yz_0_yy = buffer_pdsd[99];

    auto g_z_yz_0_yz = buffer_pdsd[100];

    auto g_z_yz_0_zz = buffer_pdsd[101];

    auto g_z_zz_0_xx = buffer_pdsd[102];

    auto g_z_zz_0_xy = buffer_pdsd[103];

    auto g_z_zz_0_xz = buffer_pdsd[104];

    auto g_z_zz_0_yy = buffer_pdsd[105];

    auto g_z_zz_0_yz = buffer_pdsd[106];

    auto g_z_zz_0_zz = buffer_pdsd[107];

    /// Set up components of auxilary buffer : buffer_fdsd

    auto g_xxx_xx_0_xx = buffer_fdsd[0];

    auto g_xxx_xx_0_xy = buffer_fdsd[1];

    auto g_xxx_xx_0_xz = buffer_fdsd[2];

    auto g_xxx_xx_0_yy = buffer_fdsd[3];

    auto g_xxx_xx_0_yz = buffer_fdsd[4];

    auto g_xxx_xx_0_zz = buffer_fdsd[5];

    auto g_xxx_xy_0_xx = buffer_fdsd[6];

    auto g_xxx_xy_0_xy = buffer_fdsd[7];

    auto g_xxx_xy_0_xz = buffer_fdsd[8];

    auto g_xxx_xy_0_yy = buffer_fdsd[9];

    auto g_xxx_xy_0_yz = buffer_fdsd[10];

    auto g_xxx_xy_0_zz = buffer_fdsd[11];

    auto g_xxx_xz_0_xx = buffer_fdsd[12];

    auto g_xxx_xz_0_xy = buffer_fdsd[13];

    auto g_xxx_xz_0_xz = buffer_fdsd[14];

    auto g_xxx_xz_0_yy = buffer_fdsd[15];

    auto g_xxx_xz_0_yz = buffer_fdsd[16];

    auto g_xxx_xz_0_zz = buffer_fdsd[17];

    auto g_xxx_yy_0_xx = buffer_fdsd[18];

    auto g_xxx_yy_0_xy = buffer_fdsd[19];

    auto g_xxx_yy_0_xz = buffer_fdsd[20];

    auto g_xxx_yy_0_yy = buffer_fdsd[21];

    auto g_xxx_yy_0_yz = buffer_fdsd[22];

    auto g_xxx_yy_0_zz = buffer_fdsd[23];

    auto g_xxx_yz_0_xx = buffer_fdsd[24];

    auto g_xxx_yz_0_xy = buffer_fdsd[25];

    auto g_xxx_yz_0_xz = buffer_fdsd[26];

    auto g_xxx_yz_0_yy = buffer_fdsd[27];

    auto g_xxx_yz_0_yz = buffer_fdsd[28];

    auto g_xxx_yz_0_zz = buffer_fdsd[29];

    auto g_xxx_zz_0_xx = buffer_fdsd[30];

    auto g_xxx_zz_0_xy = buffer_fdsd[31];

    auto g_xxx_zz_0_xz = buffer_fdsd[32];

    auto g_xxx_zz_0_yy = buffer_fdsd[33];

    auto g_xxx_zz_0_yz = buffer_fdsd[34];

    auto g_xxx_zz_0_zz = buffer_fdsd[35];

    auto g_xxy_xx_0_xx = buffer_fdsd[36];

    auto g_xxy_xx_0_xy = buffer_fdsd[37];

    auto g_xxy_xx_0_xz = buffer_fdsd[38];

    auto g_xxy_xx_0_yy = buffer_fdsd[39];

    auto g_xxy_xx_0_yz = buffer_fdsd[40];

    auto g_xxy_xx_0_zz = buffer_fdsd[41];

    auto g_xxy_xy_0_xx = buffer_fdsd[42];

    auto g_xxy_xy_0_xy = buffer_fdsd[43];

    auto g_xxy_xy_0_xz = buffer_fdsd[44];

    auto g_xxy_xy_0_yy = buffer_fdsd[45];

    auto g_xxy_xy_0_yz = buffer_fdsd[46];

    auto g_xxy_xy_0_zz = buffer_fdsd[47];

    auto g_xxy_xz_0_xx = buffer_fdsd[48];

    auto g_xxy_xz_0_xy = buffer_fdsd[49];

    auto g_xxy_xz_0_xz = buffer_fdsd[50];

    auto g_xxy_xz_0_yy = buffer_fdsd[51];

    auto g_xxy_xz_0_yz = buffer_fdsd[52];

    auto g_xxy_xz_0_zz = buffer_fdsd[53];

    auto g_xxy_yy_0_xx = buffer_fdsd[54];

    auto g_xxy_yy_0_xy = buffer_fdsd[55];

    auto g_xxy_yy_0_xz = buffer_fdsd[56];

    auto g_xxy_yy_0_yy = buffer_fdsd[57];

    auto g_xxy_yy_0_yz = buffer_fdsd[58];

    auto g_xxy_yy_0_zz = buffer_fdsd[59];

    auto g_xxy_yz_0_xx = buffer_fdsd[60];

    auto g_xxy_yz_0_xy = buffer_fdsd[61];

    auto g_xxy_yz_0_xz = buffer_fdsd[62];

    auto g_xxy_yz_0_yy = buffer_fdsd[63];

    auto g_xxy_yz_0_yz = buffer_fdsd[64];

    auto g_xxy_yz_0_zz = buffer_fdsd[65];

    auto g_xxy_zz_0_xx = buffer_fdsd[66];

    auto g_xxy_zz_0_xy = buffer_fdsd[67];

    auto g_xxy_zz_0_xz = buffer_fdsd[68];

    auto g_xxy_zz_0_yy = buffer_fdsd[69];

    auto g_xxy_zz_0_yz = buffer_fdsd[70];

    auto g_xxy_zz_0_zz = buffer_fdsd[71];

    auto g_xxz_xx_0_xx = buffer_fdsd[72];

    auto g_xxz_xx_0_xy = buffer_fdsd[73];

    auto g_xxz_xx_0_xz = buffer_fdsd[74];

    auto g_xxz_xx_0_yy = buffer_fdsd[75];

    auto g_xxz_xx_0_yz = buffer_fdsd[76];

    auto g_xxz_xx_0_zz = buffer_fdsd[77];

    auto g_xxz_xy_0_xx = buffer_fdsd[78];

    auto g_xxz_xy_0_xy = buffer_fdsd[79];

    auto g_xxz_xy_0_xz = buffer_fdsd[80];

    auto g_xxz_xy_0_yy = buffer_fdsd[81];

    auto g_xxz_xy_0_yz = buffer_fdsd[82];

    auto g_xxz_xy_0_zz = buffer_fdsd[83];

    auto g_xxz_xz_0_xx = buffer_fdsd[84];

    auto g_xxz_xz_0_xy = buffer_fdsd[85];

    auto g_xxz_xz_0_xz = buffer_fdsd[86];

    auto g_xxz_xz_0_yy = buffer_fdsd[87];

    auto g_xxz_xz_0_yz = buffer_fdsd[88];

    auto g_xxz_xz_0_zz = buffer_fdsd[89];

    auto g_xxz_yy_0_xx = buffer_fdsd[90];

    auto g_xxz_yy_0_xy = buffer_fdsd[91];

    auto g_xxz_yy_0_xz = buffer_fdsd[92];

    auto g_xxz_yy_0_yy = buffer_fdsd[93];

    auto g_xxz_yy_0_yz = buffer_fdsd[94];

    auto g_xxz_yy_0_zz = buffer_fdsd[95];

    auto g_xxz_yz_0_xx = buffer_fdsd[96];

    auto g_xxz_yz_0_xy = buffer_fdsd[97];

    auto g_xxz_yz_0_xz = buffer_fdsd[98];

    auto g_xxz_yz_0_yy = buffer_fdsd[99];

    auto g_xxz_yz_0_yz = buffer_fdsd[100];

    auto g_xxz_yz_0_zz = buffer_fdsd[101];

    auto g_xxz_zz_0_xx = buffer_fdsd[102];

    auto g_xxz_zz_0_xy = buffer_fdsd[103];

    auto g_xxz_zz_0_xz = buffer_fdsd[104];

    auto g_xxz_zz_0_yy = buffer_fdsd[105];

    auto g_xxz_zz_0_yz = buffer_fdsd[106];

    auto g_xxz_zz_0_zz = buffer_fdsd[107];

    auto g_xyy_xx_0_xx = buffer_fdsd[108];

    auto g_xyy_xx_0_xy = buffer_fdsd[109];

    auto g_xyy_xx_0_xz = buffer_fdsd[110];

    auto g_xyy_xx_0_yy = buffer_fdsd[111];

    auto g_xyy_xx_0_yz = buffer_fdsd[112];

    auto g_xyy_xx_0_zz = buffer_fdsd[113];

    auto g_xyy_xy_0_xx = buffer_fdsd[114];

    auto g_xyy_xy_0_xy = buffer_fdsd[115];

    auto g_xyy_xy_0_xz = buffer_fdsd[116];

    auto g_xyy_xy_0_yy = buffer_fdsd[117];

    auto g_xyy_xy_0_yz = buffer_fdsd[118];

    auto g_xyy_xy_0_zz = buffer_fdsd[119];

    auto g_xyy_xz_0_xx = buffer_fdsd[120];

    auto g_xyy_xz_0_xy = buffer_fdsd[121];

    auto g_xyy_xz_0_xz = buffer_fdsd[122];

    auto g_xyy_xz_0_yy = buffer_fdsd[123];

    auto g_xyy_xz_0_yz = buffer_fdsd[124];

    auto g_xyy_xz_0_zz = buffer_fdsd[125];

    auto g_xyy_yy_0_xx = buffer_fdsd[126];

    auto g_xyy_yy_0_xy = buffer_fdsd[127];

    auto g_xyy_yy_0_xz = buffer_fdsd[128];

    auto g_xyy_yy_0_yy = buffer_fdsd[129];

    auto g_xyy_yy_0_yz = buffer_fdsd[130];

    auto g_xyy_yy_0_zz = buffer_fdsd[131];

    auto g_xyy_yz_0_xx = buffer_fdsd[132];

    auto g_xyy_yz_0_xy = buffer_fdsd[133];

    auto g_xyy_yz_0_xz = buffer_fdsd[134];

    auto g_xyy_yz_0_yy = buffer_fdsd[135];

    auto g_xyy_yz_0_yz = buffer_fdsd[136];

    auto g_xyy_yz_0_zz = buffer_fdsd[137];

    auto g_xyy_zz_0_xx = buffer_fdsd[138];

    auto g_xyy_zz_0_xy = buffer_fdsd[139];

    auto g_xyy_zz_0_xz = buffer_fdsd[140];

    auto g_xyy_zz_0_yy = buffer_fdsd[141];

    auto g_xyy_zz_0_yz = buffer_fdsd[142];

    auto g_xyy_zz_0_zz = buffer_fdsd[143];

    auto g_xyz_xx_0_xx = buffer_fdsd[144];

    auto g_xyz_xx_0_xy = buffer_fdsd[145];

    auto g_xyz_xx_0_xz = buffer_fdsd[146];

    auto g_xyz_xx_0_yy = buffer_fdsd[147];

    auto g_xyz_xx_0_yz = buffer_fdsd[148];

    auto g_xyz_xx_0_zz = buffer_fdsd[149];

    auto g_xyz_xy_0_xx = buffer_fdsd[150];

    auto g_xyz_xy_0_xy = buffer_fdsd[151];

    auto g_xyz_xy_0_xz = buffer_fdsd[152];

    auto g_xyz_xy_0_yy = buffer_fdsd[153];

    auto g_xyz_xy_0_yz = buffer_fdsd[154];

    auto g_xyz_xy_0_zz = buffer_fdsd[155];

    auto g_xyz_xz_0_xx = buffer_fdsd[156];

    auto g_xyz_xz_0_xy = buffer_fdsd[157];

    auto g_xyz_xz_0_xz = buffer_fdsd[158];

    auto g_xyz_xz_0_yy = buffer_fdsd[159];

    auto g_xyz_xz_0_yz = buffer_fdsd[160];

    auto g_xyz_xz_0_zz = buffer_fdsd[161];

    auto g_xyz_yy_0_xx = buffer_fdsd[162];

    auto g_xyz_yy_0_xy = buffer_fdsd[163];

    auto g_xyz_yy_0_xz = buffer_fdsd[164];

    auto g_xyz_yy_0_yy = buffer_fdsd[165];

    auto g_xyz_yy_0_yz = buffer_fdsd[166];

    auto g_xyz_yy_0_zz = buffer_fdsd[167];

    auto g_xyz_yz_0_xx = buffer_fdsd[168];

    auto g_xyz_yz_0_xy = buffer_fdsd[169];

    auto g_xyz_yz_0_xz = buffer_fdsd[170];

    auto g_xyz_yz_0_yy = buffer_fdsd[171];

    auto g_xyz_yz_0_yz = buffer_fdsd[172];

    auto g_xyz_yz_0_zz = buffer_fdsd[173];

    auto g_xyz_zz_0_xx = buffer_fdsd[174];

    auto g_xyz_zz_0_xy = buffer_fdsd[175];

    auto g_xyz_zz_0_xz = buffer_fdsd[176];

    auto g_xyz_zz_0_yy = buffer_fdsd[177];

    auto g_xyz_zz_0_yz = buffer_fdsd[178];

    auto g_xyz_zz_0_zz = buffer_fdsd[179];

    auto g_xzz_xx_0_xx = buffer_fdsd[180];

    auto g_xzz_xx_0_xy = buffer_fdsd[181];

    auto g_xzz_xx_0_xz = buffer_fdsd[182];

    auto g_xzz_xx_0_yy = buffer_fdsd[183];

    auto g_xzz_xx_0_yz = buffer_fdsd[184];

    auto g_xzz_xx_0_zz = buffer_fdsd[185];

    auto g_xzz_xy_0_xx = buffer_fdsd[186];

    auto g_xzz_xy_0_xy = buffer_fdsd[187];

    auto g_xzz_xy_0_xz = buffer_fdsd[188];

    auto g_xzz_xy_0_yy = buffer_fdsd[189];

    auto g_xzz_xy_0_yz = buffer_fdsd[190];

    auto g_xzz_xy_0_zz = buffer_fdsd[191];

    auto g_xzz_xz_0_xx = buffer_fdsd[192];

    auto g_xzz_xz_0_xy = buffer_fdsd[193];

    auto g_xzz_xz_0_xz = buffer_fdsd[194];

    auto g_xzz_xz_0_yy = buffer_fdsd[195];

    auto g_xzz_xz_0_yz = buffer_fdsd[196];

    auto g_xzz_xz_0_zz = buffer_fdsd[197];

    auto g_xzz_yy_0_xx = buffer_fdsd[198];

    auto g_xzz_yy_0_xy = buffer_fdsd[199];

    auto g_xzz_yy_0_xz = buffer_fdsd[200];

    auto g_xzz_yy_0_yy = buffer_fdsd[201];

    auto g_xzz_yy_0_yz = buffer_fdsd[202];

    auto g_xzz_yy_0_zz = buffer_fdsd[203];

    auto g_xzz_yz_0_xx = buffer_fdsd[204];

    auto g_xzz_yz_0_xy = buffer_fdsd[205];

    auto g_xzz_yz_0_xz = buffer_fdsd[206];

    auto g_xzz_yz_0_yy = buffer_fdsd[207];

    auto g_xzz_yz_0_yz = buffer_fdsd[208];

    auto g_xzz_yz_0_zz = buffer_fdsd[209];

    auto g_xzz_zz_0_xx = buffer_fdsd[210];

    auto g_xzz_zz_0_xy = buffer_fdsd[211];

    auto g_xzz_zz_0_xz = buffer_fdsd[212];

    auto g_xzz_zz_0_yy = buffer_fdsd[213];

    auto g_xzz_zz_0_yz = buffer_fdsd[214];

    auto g_xzz_zz_0_zz = buffer_fdsd[215];

    auto g_yyy_xx_0_xx = buffer_fdsd[216];

    auto g_yyy_xx_0_xy = buffer_fdsd[217];

    auto g_yyy_xx_0_xz = buffer_fdsd[218];

    auto g_yyy_xx_0_yy = buffer_fdsd[219];

    auto g_yyy_xx_0_yz = buffer_fdsd[220];

    auto g_yyy_xx_0_zz = buffer_fdsd[221];

    auto g_yyy_xy_0_xx = buffer_fdsd[222];

    auto g_yyy_xy_0_xy = buffer_fdsd[223];

    auto g_yyy_xy_0_xz = buffer_fdsd[224];

    auto g_yyy_xy_0_yy = buffer_fdsd[225];

    auto g_yyy_xy_0_yz = buffer_fdsd[226];

    auto g_yyy_xy_0_zz = buffer_fdsd[227];

    auto g_yyy_xz_0_xx = buffer_fdsd[228];

    auto g_yyy_xz_0_xy = buffer_fdsd[229];

    auto g_yyy_xz_0_xz = buffer_fdsd[230];

    auto g_yyy_xz_0_yy = buffer_fdsd[231];

    auto g_yyy_xz_0_yz = buffer_fdsd[232];

    auto g_yyy_xz_0_zz = buffer_fdsd[233];

    auto g_yyy_yy_0_xx = buffer_fdsd[234];

    auto g_yyy_yy_0_xy = buffer_fdsd[235];

    auto g_yyy_yy_0_xz = buffer_fdsd[236];

    auto g_yyy_yy_0_yy = buffer_fdsd[237];

    auto g_yyy_yy_0_yz = buffer_fdsd[238];

    auto g_yyy_yy_0_zz = buffer_fdsd[239];

    auto g_yyy_yz_0_xx = buffer_fdsd[240];

    auto g_yyy_yz_0_xy = buffer_fdsd[241];

    auto g_yyy_yz_0_xz = buffer_fdsd[242];

    auto g_yyy_yz_0_yy = buffer_fdsd[243];

    auto g_yyy_yz_0_yz = buffer_fdsd[244];

    auto g_yyy_yz_0_zz = buffer_fdsd[245];

    auto g_yyy_zz_0_xx = buffer_fdsd[246];

    auto g_yyy_zz_0_xy = buffer_fdsd[247];

    auto g_yyy_zz_0_xz = buffer_fdsd[248];

    auto g_yyy_zz_0_yy = buffer_fdsd[249];

    auto g_yyy_zz_0_yz = buffer_fdsd[250];

    auto g_yyy_zz_0_zz = buffer_fdsd[251];

    auto g_yyz_xx_0_xx = buffer_fdsd[252];

    auto g_yyz_xx_0_xy = buffer_fdsd[253];

    auto g_yyz_xx_0_xz = buffer_fdsd[254];

    auto g_yyz_xx_0_yy = buffer_fdsd[255];

    auto g_yyz_xx_0_yz = buffer_fdsd[256];

    auto g_yyz_xx_0_zz = buffer_fdsd[257];

    auto g_yyz_xy_0_xx = buffer_fdsd[258];

    auto g_yyz_xy_0_xy = buffer_fdsd[259];

    auto g_yyz_xy_0_xz = buffer_fdsd[260];

    auto g_yyz_xy_0_yy = buffer_fdsd[261];

    auto g_yyz_xy_0_yz = buffer_fdsd[262];

    auto g_yyz_xy_0_zz = buffer_fdsd[263];

    auto g_yyz_xz_0_xx = buffer_fdsd[264];

    auto g_yyz_xz_0_xy = buffer_fdsd[265];

    auto g_yyz_xz_0_xz = buffer_fdsd[266];

    auto g_yyz_xz_0_yy = buffer_fdsd[267];

    auto g_yyz_xz_0_yz = buffer_fdsd[268];

    auto g_yyz_xz_0_zz = buffer_fdsd[269];

    auto g_yyz_yy_0_xx = buffer_fdsd[270];

    auto g_yyz_yy_0_xy = buffer_fdsd[271];

    auto g_yyz_yy_0_xz = buffer_fdsd[272];

    auto g_yyz_yy_0_yy = buffer_fdsd[273];

    auto g_yyz_yy_0_yz = buffer_fdsd[274];

    auto g_yyz_yy_0_zz = buffer_fdsd[275];

    auto g_yyz_yz_0_xx = buffer_fdsd[276];

    auto g_yyz_yz_0_xy = buffer_fdsd[277];

    auto g_yyz_yz_0_xz = buffer_fdsd[278];

    auto g_yyz_yz_0_yy = buffer_fdsd[279];

    auto g_yyz_yz_0_yz = buffer_fdsd[280];

    auto g_yyz_yz_0_zz = buffer_fdsd[281];

    auto g_yyz_zz_0_xx = buffer_fdsd[282];

    auto g_yyz_zz_0_xy = buffer_fdsd[283];

    auto g_yyz_zz_0_xz = buffer_fdsd[284];

    auto g_yyz_zz_0_yy = buffer_fdsd[285];

    auto g_yyz_zz_0_yz = buffer_fdsd[286];

    auto g_yyz_zz_0_zz = buffer_fdsd[287];

    auto g_yzz_xx_0_xx = buffer_fdsd[288];

    auto g_yzz_xx_0_xy = buffer_fdsd[289];

    auto g_yzz_xx_0_xz = buffer_fdsd[290];

    auto g_yzz_xx_0_yy = buffer_fdsd[291];

    auto g_yzz_xx_0_yz = buffer_fdsd[292];

    auto g_yzz_xx_0_zz = buffer_fdsd[293];

    auto g_yzz_xy_0_xx = buffer_fdsd[294];

    auto g_yzz_xy_0_xy = buffer_fdsd[295];

    auto g_yzz_xy_0_xz = buffer_fdsd[296];

    auto g_yzz_xy_0_yy = buffer_fdsd[297];

    auto g_yzz_xy_0_yz = buffer_fdsd[298];

    auto g_yzz_xy_0_zz = buffer_fdsd[299];

    auto g_yzz_xz_0_xx = buffer_fdsd[300];

    auto g_yzz_xz_0_xy = buffer_fdsd[301];

    auto g_yzz_xz_0_xz = buffer_fdsd[302];

    auto g_yzz_xz_0_yy = buffer_fdsd[303];

    auto g_yzz_xz_0_yz = buffer_fdsd[304];

    auto g_yzz_xz_0_zz = buffer_fdsd[305];

    auto g_yzz_yy_0_xx = buffer_fdsd[306];

    auto g_yzz_yy_0_xy = buffer_fdsd[307];

    auto g_yzz_yy_0_xz = buffer_fdsd[308];

    auto g_yzz_yy_0_yy = buffer_fdsd[309];

    auto g_yzz_yy_0_yz = buffer_fdsd[310];

    auto g_yzz_yy_0_zz = buffer_fdsd[311];

    auto g_yzz_yz_0_xx = buffer_fdsd[312];

    auto g_yzz_yz_0_xy = buffer_fdsd[313];

    auto g_yzz_yz_0_xz = buffer_fdsd[314];

    auto g_yzz_yz_0_yy = buffer_fdsd[315];

    auto g_yzz_yz_0_yz = buffer_fdsd[316];

    auto g_yzz_yz_0_zz = buffer_fdsd[317];

    auto g_yzz_zz_0_xx = buffer_fdsd[318];

    auto g_yzz_zz_0_xy = buffer_fdsd[319];

    auto g_yzz_zz_0_xz = buffer_fdsd[320];

    auto g_yzz_zz_0_yy = buffer_fdsd[321];

    auto g_yzz_zz_0_yz = buffer_fdsd[322];

    auto g_yzz_zz_0_zz = buffer_fdsd[323];

    auto g_zzz_xx_0_xx = buffer_fdsd[324];

    auto g_zzz_xx_0_xy = buffer_fdsd[325];

    auto g_zzz_xx_0_xz = buffer_fdsd[326];

    auto g_zzz_xx_0_yy = buffer_fdsd[327];

    auto g_zzz_xx_0_yz = buffer_fdsd[328];

    auto g_zzz_xx_0_zz = buffer_fdsd[329];

    auto g_zzz_xy_0_xx = buffer_fdsd[330];

    auto g_zzz_xy_0_xy = buffer_fdsd[331];

    auto g_zzz_xy_0_xz = buffer_fdsd[332];

    auto g_zzz_xy_0_yy = buffer_fdsd[333];

    auto g_zzz_xy_0_yz = buffer_fdsd[334];

    auto g_zzz_xy_0_zz = buffer_fdsd[335];

    auto g_zzz_xz_0_xx = buffer_fdsd[336];

    auto g_zzz_xz_0_xy = buffer_fdsd[337];

    auto g_zzz_xz_0_xz = buffer_fdsd[338];

    auto g_zzz_xz_0_yy = buffer_fdsd[339];

    auto g_zzz_xz_0_yz = buffer_fdsd[340];

    auto g_zzz_xz_0_zz = buffer_fdsd[341];

    auto g_zzz_yy_0_xx = buffer_fdsd[342];

    auto g_zzz_yy_0_xy = buffer_fdsd[343];

    auto g_zzz_yy_0_xz = buffer_fdsd[344];

    auto g_zzz_yy_0_yy = buffer_fdsd[345];

    auto g_zzz_yy_0_yz = buffer_fdsd[346];

    auto g_zzz_yy_0_zz = buffer_fdsd[347];

    auto g_zzz_yz_0_xx = buffer_fdsd[348];

    auto g_zzz_yz_0_xy = buffer_fdsd[349];

    auto g_zzz_yz_0_xz = buffer_fdsd[350];

    auto g_zzz_yz_0_yy = buffer_fdsd[351];

    auto g_zzz_yz_0_yz = buffer_fdsd[352];

    auto g_zzz_yz_0_zz = buffer_fdsd[353];

    auto g_zzz_zz_0_xx = buffer_fdsd[354];

    auto g_zzz_zz_0_xy = buffer_fdsd[355];

    auto g_zzz_zz_0_xz = buffer_fdsd[356];

    auto g_zzz_zz_0_yy = buffer_fdsd[357];

    auto g_zzz_zz_0_yz = buffer_fdsd[358];

    auto g_zzz_zz_0_zz = buffer_fdsd[359];

    /// Set up components of integrals buffer : buffer_2000_pdsd

    auto g_xx_0_0_0_x_xx_0_xx = buffer_2000_pdsd[0];

    auto g_xx_0_0_0_x_xx_0_xy = buffer_2000_pdsd[1];

    auto g_xx_0_0_0_x_xx_0_xz = buffer_2000_pdsd[2];

    auto g_xx_0_0_0_x_xx_0_yy = buffer_2000_pdsd[3];

    auto g_xx_0_0_0_x_xx_0_yz = buffer_2000_pdsd[4];

    auto g_xx_0_0_0_x_xx_0_zz = buffer_2000_pdsd[5];

    auto g_xx_0_0_0_x_xy_0_xx = buffer_2000_pdsd[6];

    auto g_xx_0_0_0_x_xy_0_xy = buffer_2000_pdsd[7];

    auto g_xx_0_0_0_x_xy_0_xz = buffer_2000_pdsd[8];

    auto g_xx_0_0_0_x_xy_0_yy = buffer_2000_pdsd[9];

    auto g_xx_0_0_0_x_xy_0_yz = buffer_2000_pdsd[10];

    auto g_xx_0_0_0_x_xy_0_zz = buffer_2000_pdsd[11];

    auto g_xx_0_0_0_x_xz_0_xx = buffer_2000_pdsd[12];

    auto g_xx_0_0_0_x_xz_0_xy = buffer_2000_pdsd[13];

    auto g_xx_0_0_0_x_xz_0_xz = buffer_2000_pdsd[14];

    auto g_xx_0_0_0_x_xz_0_yy = buffer_2000_pdsd[15];

    auto g_xx_0_0_0_x_xz_0_yz = buffer_2000_pdsd[16];

    auto g_xx_0_0_0_x_xz_0_zz = buffer_2000_pdsd[17];

    auto g_xx_0_0_0_x_yy_0_xx = buffer_2000_pdsd[18];

    auto g_xx_0_0_0_x_yy_0_xy = buffer_2000_pdsd[19];

    auto g_xx_0_0_0_x_yy_0_xz = buffer_2000_pdsd[20];

    auto g_xx_0_0_0_x_yy_0_yy = buffer_2000_pdsd[21];

    auto g_xx_0_0_0_x_yy_0_yz = buffer_2000_pdsd[22];

    auto g_xx_0_0_0_x_yy_0_zz = buffer_2000_pdsd[23];

    auto g_xx_0_0_0_x_yz_0_xx = buffer_2000_pdsd[24];

    auto g_xx_0_0_0_x_yz_0_xy = buffer_2000_pdsd[25];

    auto g_xx_0_0_0_x_yz_0_xz = buffer_2000_pdsd[26];

    auto g_xx_0_0_0_x_yz_0_yy = buffer_2000_pdsd[27];

    auto g_xx_0_0_0_x_yz_0_yz = buffer_2000_pdsd[28];

    auto g_xx_0_0_0_x_yz_0_zz = buffer_2000_pdsd[29];

    auto g_xx_0_0_0_x_zz_0_xx = buffer_2000_pdsd[30];

    auto g_xx_0_0_0_x_zz_0_xy = buffer_2000_pdsd[31];

    auto g_xx_0_0_0_x_zz_0_xz = buffer_2000_pdsd[32];

    auto g_xx_0_0_0_x_zz_0_yy = buffer_2000_pdsd[33];

    auto g_xx_0_0_0_x_zz_0_yz = buffer_2000_pdsd[34];

    auto g_xx_0_0_0_x_zz_0_zz = buffer_2000_pdsd[35];

    auto g_xx_0_0_0_y_xx_0_xx = buffer_2000_pdsd[36];

    auto g_xx_0_0_0_y_xx_0_xy = buffer_2000_pdsd[37];

    auto g_xx_0_0_0_y_xx_0_xz = buffer_2000_pdsd[38];

    auto g_xx_0_0_0_y_xx_0_yy = buffer_2000_pdsd[39];

    auto g_xx_0_0_0_y_xx_0_yz = buffer_2000_pdsd[40];

    auto g_xx_0_0_0_y_xx_0_zz = buffer_2000_pdsd[41];

    auto g_xx_0_0_0_y_xy_0_xx = buffer_2000_pdsd[42];

    auto g_xx_0_0_0_y_xy_0_xy = buffer_2000_pdsd[43];

    auto g_xx_0_0_0_y_xy_0_xz = buffer_2000_pdsd[44];

    auto g_xx_0_0_0_y_xy_0_yy = buffer_2000_pdsd[45];

    auto g_xx_0_0_0_y_xy_0_yz = buffer_2000_pdsd[46];

    auto g_xx_0_0_0_y_xy_0_zz = buffer_2000_pdsd[47];

    auto g_xx_0_0_0_y_xz_0_xx = buffer_2000_pdsd[48];

    auto g_xx_0_0_0_y_xz_0_xy = buffer_2000_pdsd[49];

    auto g_xx_0_0_0_y_xz_0_xz = buffer_2000_pdsd[50];

    auto g_xx_0_0_0_y_xz_0_yy = buffer_2000_pdsd[51];

    auto g_xx_0_0_0_y_xz_0_yz = buffer_2000_pdsd[52];

    auto g_xx_0_0_0_y_xz_0_zz = buffer_2000_pdsd[53];

    auto g_xx_0_0_0_y_yy_0_xx = buffer_2000_pdsd[54];

    auto g_xx_0_0_0_y_yy_0_xy = buffer_2000_pdsd[55];

    auto g_xx_0_0_0_y_yy_0_xz = buffer_2000_pdsd[56];

    auto g_xx_0_0_0_y_yy_0_yy = buffer_2000_pdsd[57];

    auto g_xx_0_0_0_y_yy_0_yz = buffer_2000_pdsd[58];

    auto g_xx_0_0_0_y_yy_0_zz = buffer_2000_pdsd[59];

    auto g_xx_0_0_0_y_yz_0_xx = buffer_2000_pdsd[60];

    auto g_xx_0_0_0_y_yz_0_xy = buffer_2000_pdsd[61];

    auto g_xx_0_0_0_y_yz_0_xz = buffer_2000_pdsd[62];

    auto g_xx_0_0_0_y_yz_0_yy = buffer_2000_pdsd[63];

    auto g_xx_0_0_0_y_yz_0_yz = buffer_2000_pdsd[64];

    auto g_xx_0_0_0_y_yz_0_zz = buffer_2000_pdsd[65];

    auto g_xx_0_0_0_y_zz_0_xx = buffer_2000_pdsd[66];

    auto g_xx_0_0_0_y_zz_0_xy = buffer_2000_pdsd[67];

    auto g_xx_0_0_0_y_zz_0_xz = buffer_2000_pdsd[68];

    auto g_xx_0_0_0_y_zz_0_yy = buffer_2000_pdsd[69];

    auto g_xx_0_0_0_y_zz_0_yz = buffer_2000_pdsd[70];

    auto g_xx_0_0_0_y_zz_0_zz = buffer_2000_pdsd[71];

    auto g_xx_0_0_0_z_xx_0_xx = buffer_2000_pdsd[72];

    auto g_xx_0_0_0_z_xx_0_xy = buffer_2000_pdsd[73];

    auto g_xx_0_0_0_z_xx_0_xz = buffer_2000_pdsd[74];

    auto g_xx_0_0_0_z_xx_0_yy = buffer_2000_pdsd[75];

    auto g_xx_0_0_0_z_xx_0_yz = buffer_2000_pdsd[76];

    auto g_xx_0_0_0_z_xx_0_zz = buffer_2000_pdsd[77];

    auto g_xx_0_0_0_z_xy_0_xx = buffer_2000_pdsd[78];

    auto g_xx_0_0_0_z_xy_0_xy = buffer_2000_pdsd[79];

    auto g_xx_0_0_0_z_xy_0_xz = buffer_2000_pdsd[80];

    auto g_xx_0_0_0_z_xy_0_yy = buffer_2000_pdsd[81];

    auto g_xx_0_0_0_z_xy_0_yz = buffer_2000_pdsd[82];

    auto g_xx_0_0_0_z_xy_0_zz = buffer_2000_pdsd[83];

    auto g_xx_0_0_0_z_xz_0_xx = buffer_2000_pdsd[84];

    auto g_xx_0_0_0_z_xz_0_xy = buffer_2000_pdsd[85];

    auto g_xx_0_0_0_z_xz_0_xz = buffer_2000_pdsd[86];

    auto g_xx_0_0_0_z_xz_0_yy = buffer_2000_pdsd[87];

    auto g_xx_0_0_0_z_xz_0_yz = buffer_2000_pdsd[88];

    auto g_xx_0_0_0_z_xz_0_zz = buffer_2000_pdsd[89];

    auto g_xx_0_0_0_z_yy_0_xx = buffer_2000_pdsd[90];

    auto g_xx_0_0_0_z_yy_0_xy = buffer_2000_pdsd[91];

    auto g_xx_0_0_0_z_yy_0_xz = buffer_2000_pdsd[92];

    auto g_xx_0_0_0_z_yy_0_yy = buffer_2000_pdsd[93];

    auto g_xx_0_0_0_z_yy_0_yz = buffer_2000_pdsd[94];

    auto g_xx_0_0_0_z_yy_0_zz = buffer_2000_pdsd[95];

    auto g_xx_0_0_0_z_yz_0_xx = buffer_2000_pdsd[96];

    auto g_xx_0_0_0_z_yz_0_xy = buffer_2000_pdsd[97];

    auto g_xx_0_0_0_z_yz_0_xz = buffer_2000_pdsd[98];

    auto g_xx_0_0_0_z_yz_0_yy = buffer_2000_pdsd[99];

    auto g_xx_0_0_0_z_yz_0_yz = buffer_2000_pdsd[100];

    auto g_xx_0_0_0_z_yz_0_zz = buffer_2000_pdsd[101];

    auto g_xx_0_0_0_z_zz_0_xx = buffer_2000_pdsd[102];

    auto g_xx_0_0_0_z_zz_0_xy = buffer_2000_pdsd[103];

    auto g_xx_0_0_0_z_zz_0_xz = buffer_2000_pdsd[104];

    auto g_xx_0_0_0_z_zz_0_yy = buffer_2000_pdsd[105];

    auto g_xx_0_0_0_z_zz_0_yz = buffer_2000_pdsd[106];

    auto g_xx_0_0_0_z_zz_0_zz = buffer_2000_pdsd[107];

    auto g_xy_0_0_0_x_xx_0_xx = buffer_2000_pdsd[108];

    auto g_xy_0_0_0_x_xx_0_xy = buffer_2000_pdsd[109];

    auto g_xy_0_0_0_x_xx_0_xz = buffer_2000_pdsd[110];

    auto g_xy_0_0_0_x_xx_0_yy = buffer_2000_pdsd[111];

    auto g_xy_0_0_0_x_xx_0_yz = buffer_2000_pdsd[112];

    auto g_xy_0_0_0_x_xx_0_zz = buffer_2000_pdsd[113];

    auto g_xy_0_0_0_x_xy_0_xx = buffer_2000_pdsd[114];

    auto g_xy_0_0_0_x_xy_0_xy = buffer_2000_pdsd[115];

    auto g_xy_0_0_0_x_xy_0_xz = buffer_2000_pdsd[116];

    auto g_xy_0_0_0_x_xy_0_yy = buffer_2000_pdsd[117];

    auto g_xy_0_0_0_x_xy_0_yz = buffer_2000_pdsd[118];

    auto g_xy_0_0_0_x_xy_0_zz = buffer_2000_pdsd[119];

    auto g_xy_0_0_0_x_xz_0_xx = buffer_2000_pdsd[120];

    auto g_xy_0_0_0_x_xz_0_xy = buffer_2000_pdsd[121];

    auto g_xy_0_0_0_x_xz_0_xz = buffer_2000_pdsd[122];

    auto g_xy_0_0_0_x_xz_0_yy = buffer_2000_pdsd[123];

    auto g_xy_0_0_0_x_xz_0_yz = buffer_2000_pdsd[124];

    auto g_xy_0_0_0_x_xz_0_zz = buffer_2000_pdsd[125];

    auto g_xy_0_0_0_x_yy_0_xx = buffer_2000_pdsd[126];

    auto g_xy_0_0_0_x_yy_0_xy = buffer_2000_pdsd[127];

    auto g_xy_0_0_0_x_yy_0_xz = buffer_2000_pdsd[128];

    auto g_xy_0_0_0_x_yy_0_yy = buffer_2000_pdsd[129];

    auto g_xy_0_0_0_x_yy_0_yz = buffer_2000_pdsd[130];

    auto g_xy_0_0_0_x_yy_0_zz = buffer_2000_pdsd[131];

    auto g_xy_0_0_0_x_yz_0_xx = buffer_2000_pdsd[132];

    auto g_xy_0_0_0_x_yz_0_xy = buffer_2000_pdsd[133];

    auto g_xy_0_0_0_x_yz_0_xz = buffer_2000_pdsd[134];

    auto g_xy_0_0_0_x_yz_0_yy = buffer_2000_pdsd[135];

    auto g_xy_0_0_0_x_yz_0_yz = buffer_2000_pdsd[136];

    auto g_xy_0_0_0_x_yz_0_zz = buffer_2000_pdsd[137];

    auto g_xy_0_0_0_x_zz_0_xx = buffer_2000_pdsd[138];

    auto g_xy_0_0_0_x_zz_0_xy = buffer_2000_pdsd[139];

    auto g_xy_0_0_0_x_zz_0_xz = buffer_2000_pdsd[140];

    auto g_xy_0_0_0_x_zz_0_yy = buffer_2000_pdsd[141];

    auto g_xy_0_0_0_x_zz_0_yz = buffer_2000_pdsd[142];

    auto g_xy_0_0_0_x_zz_0_zz = buffer_2000_pdsd[143];

    auto g_xy_0_0_0_y_xx_0_xx = buffer_2000_pdsd[144];

    auto g_xy_0_0_0_y_xx_0_xy = buffer_2000_pdsd[145];

    auto g_xy_0_0_0_y_xx_0_xz = buffer_2000_pdsd[146];

    auto g_xy_0_0_0_y_xx_0_yy = buffer_2000_pdsd[147];

    auto g_xy_0_0_0_y_xx_0_yz = buffer_2000_pdsd[148];

    auto g_xy_0_0_0_y_xx_0_zz = buffer_2000_pdsd[149];

    auto g_xy_0_0_0_y_xy_0_xx = buffer_2000_pdsd[150];

    auto g_xy_0_0_0_y_xy_0_xy = buffer_2000_pdsd[151];

    auto g_xy_0_0_0_y_xy_0_xz = buffer_2000_pdsd[152];

    auto g_xy_0_0_0_y_xy_0_yy = buffer_2000_pdsd[153];

    auto g_xy_0_0_0_y_xy_0_yz = buffer_2000_pdsd[154];

    auto g_xy_0_0_0_y_xy_0_zz = buffer_2000_pdsd[155];

    auto g_xy_0_0_0_y_xz_0_xx = buffer_2000_pdsd[156];

    auto g_xy_0_0_0_y_xz_0_xy = buffer_2000_pdsd[157];

    auto g_xy_0_0_0_y_xz_0_xz = buffer_2000_pdsd[158];

    auto g_xy_0_0_0_y_xz_0_yy = buffer_2000_pdsd[159];

    auto g_xy_0_0_0_y_xz_0_yz = buffer_2000_pdsd[160];

    auto g_xy_0_0_0_y_xz_0_zz = buffer_2000_pdsd[161];

    auto g_xy_0_0_0_y_yy_0_xx = buffer_2000_pdsd[162];

    auto g_xy_0_0_0_y_yy_0_xy = buffer_2000_pdsd[163];

    auto g_xy_0_0_0_y_yy_0_xz = buffer_2000_pdsd[164];

    auto g_xy_0_0_0_y_yy_0_yy = buffer_2000_pdsd[165];

    auto g_xy_0_0_0_y_yy_0_yz = buffer_2000_pdsd[166];

    auto g_xy_0_0_0_y_yy_0_zz = buffer_2000_pdsd[167];

    auto g_xy_0_0_0_y_yz_0_xx = buffer_2000_pdsd[168];

    auto g_xy_0_0_0_y_yz_0_xy = buffer_2000_pdsd[169];

    auto g_xy_0_0_0_y_yz_0_xz = buffer_2000_pdsd[170];

    auto g_xy_0_0_0_y_yz_0_yy = buffer_2000_pdsd[171];

    auto g_xy_0_0_0_y_yz_0_yz = buffer_2000_pdsd[172];

    auto g_xy_0_0_0_y_yz_0_zz = buffer_2000_pdsd[173];

    auto g_xy_0_0_0_y_zz_0_xx = buffer_2000_pdsd[174];

    auto g_xy_0_0_0_y_zz_0_xy = buffer_2000_pdsd[175];

    auto g_xy_0_0_0_y_zz_0_xz = buffer_2000_pdsd[176];

    auto g_xy_0_0_0_y_zz_0_yy = buffer_2000_pdsd[177];

    auto g_xy_0_0_0_y_zz_0_yz = buffer_2000_pdsd[178];

    auto g_xy_0_0_0_y_zz_0_zz = buffer_2000_pdsd[179];

    auto g_xy_0_0_0_z_xx_0_xx = buffer_2000_pdsd[180];

    auto g_xy_0_0_0_z_xx_0_xy = buffer_2000_pdsd[181];

    auto g_xy_0_0_0_z_xx_0_xz = buffer_2000_pdsd[182];

    auto g_xy_0_0_0_z_xx_0_yy = buffer_2000_pdsd[183];

    auto g_xy_0_0_0_z_xx_0_yz = buffer_2000_pdsd[184];

    auto g_xy_0_0_0_z_xx_0_zz = buffer_2000_pdsd[185];

    auto g_xy_0_0_0_z_xy_0_xx = buffer_2000_pdsd[186];

    auto g_xy_0_0_0_z_xy_0_xy = buffer_2000_pdsd[187];

    auto g_xy_0_0_0_z_xy_0_xz = buffer_2000_pdsd[188];

    auto g_xy_0_0_0_z_xy_0_yy = buffer_2000_pdsd[189];

    auto g_xy_0_0_0_z_xy_0_yz = buffer_2000_pdsd[190];

    auto g_xy_0_0_0_z_xy_0_zz = buffer_2000_pdsd[191];

    auto g_xy_0_0_0_z_xz_0_xx = buffer_2000_pdsd[192];

    auto g_xy_0_0_0_z_xz_0_xy = buffer_2000_pdsd[193];

    auto g_xy_0_0_0_z_xz_0_xz = buffer_2000_pdsd[194];

    auto g_xy_0_0_0_z_xz_0_yy = buffer_2000_pdsd[195];

    auto g_xy_0_0_0_z_xz_0_yz = buffer_2000_pdsd[196];

    auto g_xy_0_0_0_z_xz_0_zz = buffer_2000_pdsd[197];

    auto g_xy_0_0_0_z_yy_0_xx = buffer_2000_pdsd[198];

    auto g_xy_0_0_0_z_yy_0_xy = buffer_2000_pdsd[199];

    auto g_xy_0_0_0_z_yy_0_xz = buffer_2000_pdsd[200];

    auto g_xy_0_0_0_z_yy_0_yy = buffer_2000_pdsd[201];

    auto g_xy_0_0_0_z_yy_0_yz = buffer_2000_pdsd[202];

    auto g_xy_0_0_0_z_yy_0_zz = buffer_2000_pdsd[203];

    auto g_xy_0_0_0_z_yz_0_xx = buffer_2000_pdsd[204];

    auto g_xy_0_0_0_z_yz_0_xy = buffer_2000_pdsd[205];

    auto g_xy_0_0_0_z_yz_0_xz = buffer_2000_pdsd[206];

    auto g_xy_0_0_0_z_yz_0_yy = buffer_2000_pdsd[207];

    auto g_xy_0_0_0_z_yz_0_yz = buffer_2000_pdsd[208];

    auto g_xy_0_0_0_z_yz_0_zz = buffer_2000_pdsd[209];

    auto g_xy_0_0_0_z_zz_0_xx = buffer_2000_pdsd[210];

    auto g_xy_0_0_0_z_zz_0_xy = buffer_2000_pdsd[211];

    auto g_xy_0_0_0_z_zz_0_xz = buffer_2000_pdsd[212];

    auto g_xy_0_0_0_z_zz_0_yy = buffer_2000_pdsd[213];

    auto g_xy_0_0_0_z_zz_0_yz = buffer_2000_pdsd[214];

    auto g_xy_0_0_0_z_zz_0_zz = buffer_2000_pdsd[215];

    auto g_xz_0_0_0_x_xx_0_xx = buffer_2000_pdsd[216];

    auto g_xz_0_0_0_x_xx_0_xy = buffer_2000_pdsd[217];

    auto g_xz_0_0_0_x_xx_0_xz = buffer_2000_pdsd[218];

    auto g_xz_0_0_0_x_xx_0_yy = buffer_2000_pdsd[219];

    auto g_xz_0_0_0_x_xx_0_yz = buffer_2000_pdsd[220];

    auto g_xz_0_0_0_x_xx_0_zz = buffer_2000_pdsd[221];

    auto g_xz_0_0_0_x_xy_0_xx = buffer_2000_pdsd[222];

    auto g_xz_0_0_0_x_xy_0_xy = buffer_2000_pdsd[223];

    auto g_xz_0_0_0_x_xy_0_xz = buffer_2000_pdsd[224];

    auto g_xz_0_0_0_x_xy_0_yy = buffer_2000_pdsd[225];

    auto g_xz_0_0_0_x_xy_0_yz = buffer_2000_pdsd[226];

    auto g_xz_0_0_0_x_xy_0_zz = buffer_2000_pdsd[227];

    auto g_xz_0_0_0_x_xz_0_xx = buffer_2000_pdsd[228];

    auto g_xz_0_0_0_x_xz_0_xy = buffer_2000_pdsd[229];

    auto g_xz_0_0_0_x_xz_0_xz = buffer_2000_pdsd[230];

    auto g_xz_0_0_0_x_xz_0_yy = buffer_2000_pdsd[231];

    auto g_xz_0_0_0_x_xz_0_yz = buffer_2000_pdsd[232];

    auto g_xz_0_0_0_x_xz_0_zz = buffer_2000_pdsd[233];

    auto g_xz_0_0_0_x_yy_0_xx = buffer_2000_pdsd[234];

    auto g_xz_0_0_0_x_yy_0_xy = buffer_2000_pdsd[235];

    auto g_xz_0_0_0_x_yy_0_xz = buffer_2000_pdsd[236];

    auto g_xz_0_0_0_x_yy_0_yy = buffer_2000_pdsd[237];

    auto g_xz_0_0_0_x_yy_0_yz = buffer_2000_pdsd[238];

    auto g_xz_0_0_0_x_yy_0_zz = buffer_2000_pdsd[239];

    auto g_xz_0_0_0_x_yz_0_xx = buffer_2000_pdsd[240];

    auto g_xz_0_0_0_x_yz_0_xy = buffer_2000_pdsd[241];

    auto g_xz_0_0_0_x_yz_0_xz = buffer_2000_pdsd[242];

    auto g_xz_0_0_0_x_yz_0_yy = buffer_2000_pdsd[243];

    auto g_xz_0_0_0_x_yz_0_yz = buffer_2000_pdsd[244];

    auto g_xz_0_0_0_x_yz_0_zz = buffer_2000_pdsd[245];

    auto g_xz_0_0_0_x_zz_0_xx = buffer_2000_pdsd[246];

    auto g_xz_0_0_0_x_zz_0_xy = buffer_2000_pdsd[247];

    auto g_xz_0_0_0_x_zz_0_xz = buffer_2000_pdsd[248];

    auto g_xz_0_0_0_x_zz_0_yy = buffer_2000_pdsd[249];

    auto g_xz_0_0_0_x_zz_0_yz = buffer_2000_pdsd[250];

    auto g_xz_0_0_0_x_zz_0_zz = buffer_2000_pdsd[251];

    auto g_xz_0_0_0_y_xx_0_xx = buffer_2000_pdsd[252];

    auto g_xz_0_0_0_y_xx_0_xy = buffer_2000_pdsd[253];

    auto g_xz_0_0_0_y_xx_0_xz = buffer_2000_pdsd[254];

    auto g_xz_0_0_0_y_xx_0_yy = buffer_2000_pdsd[255];

    auto g_xz_0_0_0_y_xx_0_yz = buffer_2000_pdsd[256];

    auto g_xz_0_0_0_y_xx_0_zz = buffer_2000_pdsd[257];

    auto g_xz_0_0_0_y_xy_0_xx = buffer_2000_pdsd[258];

    auto g_xz_0_0_0_y_xy_0_xy = buffer_2000_pdsd[259];

    auto g_xz_0_0_0_y_xy_0_xz = buffer_2000_pdsd[260];

    auto g_xz_0_0_0_y_xy_0_yy = buffer_2000_pdsd[261];

    auto g_xz_0_0_0_y_xy_0_yz = buffer_2000_pdsd[262];

    auto g_xz_0_0_0_y_xy_0_zz = buffer_2000_pdsd[263];

    auto g_xz_0_0_0_y_xz_0_xx = buffer_2000_pdsd[264];

    auto g_xz_0_0_0_y_xz_0_xy = buffer_2000_pdsd[265];

    auto g_xz_0_0_0_y_xz_0_xz = buffer_2000_pdsd[266];

    auto g_xz_0_0_0_y_xz_0_yy = buffer_2000_pdsd[267];

    auto g_xz_0_0_0_y_xz_0_yz = buffer_2000_pdsd[268];

    auto g_xz_0_0_0_y_xz_0_zz = buffer_2000_pdsd[269];

    auto g_xz_0_0_0_y_yy_0_xx = buffer_2000_pdsd[270];

    auto g_xz_0_0_0_y_yy_0_xy = buffer_2000_pdsd[271];

    auto g_xz_0_0_0_y_yy_0_xz = buffer_2000_pdsd[272];

    auto g_xz_0_0_0_y_yy_0_yy = buffer_2000_pdsd[273];

    auto g_xz_0_0_0_y_yy_0_yz = buffer_2000_pdsd[274];

    auto g_xz_0_0_0_y_yy_0_zz = buffer_2000_pdsd[275];

    auto g_xz_0_0_0_y_yz_0_xx = buffer_2000_pdsd[276];

    auto g_xz_0_0_0_y_yz_0_xy = buffer_2000_pdsd[277];

    auto g_xz_0_0_0_y_yz_0_xz = buffer_2000_pdsd[278];

    auto g_xz_0_0_0_y_yz_0_yy = buffer_2000_pdsd[279];

    auto g_xz_0_0_0_y_yz_0_yz = buffer_2000_pdsd[280];

    auto g_xz_0_0_0_y_yz_0_zz = buffer_2000_pdsd[281];

    auto g_xz_0_0_0_y_zz_0_xx = buffer_2000_pdsd[282];

    auto g_xz_0_0_0_y_zz_0_xy = buffer_2000_pdsd[283];

    auto g_xz_0_0_0_y_zz_0_xz = buffer_2000_pdsd[284];

    auto g_xz_0_0_0_y_zz_0_yy = buffer_2000_pdsd[285];

    auto g_xz_0_0_0_y_zz_0_yz = buffer_2000_pdsd[286];

    auto g_xz_0_0_0_y_zz_0_zz = buffer_2000_pdsd[287];

    auto g_xz_0_0_0_z_xx_0_xx = buffer_2000_pdsd[288];

    auto g_xz_0_0_0_z_xx_0_xy = buffer_2000_pdsd[289];

    auto g_xz_0_0_0_z_xx_0_xz = buffer_2000_pdsd[290];

    auto g_xz_0_0_0_z_xx_0_yy = buffer_2000_pdsd[291];

    auto g_xz_0_0_0_z_xx_0_yz = buffer_2000_pdsd[292];

    auto g_xz_0_0_0_z_xx_0_zz = buffer_2000_pdsd[293];

    auto g_xz_0_0_0_z_xy_0_xx = buffer_2000_pdsd[294];

    auto g_xz_0_0_0_z_xy_0_xy = buffer_2000_pdsd[295];

    auto g_xz_0_0_0_z_xy_0_xz = buffer_2000_pdsd[296];

    auto g_xz_0_0_0_z_xy_0_yy = buffer_2000_pdsd[297];

    auto g_xz_0_0_0_z_xy_0_yz = buffer_2000_pdsd[298];

    auto g_xz_0_0_0_z_xy_0_zz = buffer_2000_pdsd[299];

    auto g_xz_0_0_0_z_xz_0_xx = buffer_2000_pdsd[300];

    auto g_xz_0_0_0_z_xz_0_xy = buffer_2000_pdsd[301];

    auto g_xz_0_0_0_z_xz_0_xz = buffer_2000_pdsd[302];

    auto g_xz_0_0_0_z_xz_0_yy = buffer_2000_pdsd[303];

    auto g_xz_0_0_0_z_xz_0_yz = buffer_2000_pdsd[304];

    auto g_xz_0_0_0_z_xz_0_zz = buffer_2000_pdsd[305];

    auto g_xz_0_0_0_z_yy_0_xx = buffer_2000_pdsd[306];

    auto g_xz_0_0_0_z_yy_0_xy = buffer_2000_pdsd[307];

    auto g_xz_0_0_0_z_yy_0_xz = buffer_2000_pdsd[308];

    auto g_xz_0_0_0_z_yy_0_yy = buffer_2000_pdsd[309];

    auto g_xz_0_0_0_z_yy_0_yz = buffer_2000_pdsd[310];

    auto g_xz_0_0_0_z_yy_0_zz = buffer_2000_pdsd[311];

    auto g_xz_0_0_0_z_yz_0_xx = buffer_2000_pdsd[312];

    auto g_xz_0_0_0_z_yz_0_xy = buffer_2000_pdsd[313];

    auto g_xz_0_0_0_z_yz_0_xz = buffer_2000_pdsd[314];

    auto g_xz_0_0_0_z_yz_0_yy = buffer_2000_pdsd[315];

    auto g_xz_0_0_0_z_yz_0_yz = buffer_2000_pdsd[316];

    auto g_xz_0_0_0_z_yz_0_zz = buffer_2000_pdsd[317];

    auto g_xz_0_0_0_z_zz_0_xx = buffer_2000_pdsd[318];

    auto g_xz_0_0_0_z_zz_0_xy = buffer_2000_pdsd[319];

    auto g_xz_0_0_0_z_zz_0_xz = buffer_2000_pdsd[320];

    auto g_xz_0_0_0_z_zz_0_yy = buffer_2000_pdsd[321];

    auto g_xz_0_0_0_z_zz_0_yz = buffer_2000_pdsd[322];

    auto g_xz_0_0_0_z_zz_0_zz = buffer_2000_pdsd[323];

    auto g_yy_0_0_0_x_xx_0_xx = buffer_2000_pdsd[324];

    auto g_yy_0_0_0_x_xx_0_xy = buffer_2000_pdsd[325];

    auto g_yy_0_0_0_x_xx_0_xz = buffer_2000_pdsd[326];

    auto g_yy_0_0_0_x_xx_0_yy = buffer_2000_pdsd[327];

    auto g_yy_0_0_0_x_xx_0_yz = buffer_2000_pdsd[328];

    auto g_yy_0_0_0_x_xx_0_zz = buffer_2000_pdsd[329];

    auto g_yy_0_0_0_x_xy_0_xx = buffer_2000_pdsd[330];

    auto g_yy_0_0_0_x_xy_0_xy = buffer_2000_pdsd[331];

    auto g_yy_0_0_0_x_xy_0_xz = buffer_2000_pdsd[332];

    auto g_yy_0_0_0_x_xy_0_yy = buffer_2000_pdsd[333];

    auto g_yy_0_0_0_x_xy_0_yz = buffer_2000_pdsd[334];

    auto g_yy_0_0_0_x_xy_0_zz = buffer_2000_pdsd[335];

    auto g_yy_0_0_0_x_xz_0_xx = buffer_2000_pdsd[336];

    auto g_yy_0_0_0_x_xz_0_xy = buffer_2000_pdsd[337];

    auto g_yy_0_0_0_x_xz_0_xz = buffer_2000_pdsd[338];

    auto g_yy_0_0_0_x_xz_0_yy = buffer_2000_pdsd[339];

    auto g_yy_0_0_0_x_xz_0_yz = buffer_2000_pdsd[340];

    auto g_yy_0_0_0_x_xz_0_zz = buffer_2000_pdsd[341];

    auto g_yy_0_0_0_x_yy_0_xx = buffer_2000_pdsd[342];

    auto g_yy_0_0_0_x_yy_0_xy = buffer_2000_pdsd[343];

    auto g_yy_0_0_0_x_yy_0_xz = buffer_2000_pdsd[344];

    auto g_yy_0_0_0_x_yy_0_yy = buffer_2000_pdsd[345];

    auto g_yy_0_0_0_x_yy_0_yz = buffer_2000_pdsd[346];

    auto g_yy_0_0_0_x_yy_0_zz = buffer_2000_pdsd[347];

    auto g_yy_0_0_0_x_yz_0_xx = buffer_2000_pdsd[348];

    auto g_yy_0_0_0_x_yz_0_xy = buffer_2000_pdsd[349];

    auto g_yy_0_0_0_x_yz_0_xz = buffer_2000_pdsd[350];

    auto g_yy_0_0_0_x_yz_0_yy = buffer_2000_pdsd[351];

    auto g_yy_0_0_0_x_yz_0_yz = buffer_2000_pdsd[352];

    auto g_yy_0_0_0_x_yz_0_zz = buffer_2000_pdsd[353];

    auto g_yy_0_0_0_x_zz_0_xx = buffer_2000_pdsd[354];

    auto g_yy_0_0_0_x_zz_0_xy = buffer_2000_pdsd[355];

    auto g_yy_0_0_0_x_zz_0_xz = buffer_2000_pdsd[356];

    auto g_yy_0_0_0_x_zz_0_yy = buffer_2000_pdsd[357];

    auto g_yy_0_0_0_x_zz_0_yz = buffer_2000_pdsd[358];

    auto g_yy_0_0_0_x_zz_0_zz = buffer_2000_pdsd[359];

    auto g_yy_0_0_0_y_xx_0_xx = buffer_2000_pdsd[360];

    auto g_yy_0_0_0_y_xx_0_xy = buffer_2000_pdsd[361];

    auto g_yy_0_0_0_y_xx_0_xz = buffer_2000_pdsd[362];

    auto g_yy_0_0_0_y_xx_0_yy = buffer_2000_pdsd[363];

    auto g_yy_0_0_0_y_xx_0_yz = buffer_2000_pdsd[364];

    auto g_yy_0_0_0_y_xx_0_zz = buffer_2000_pdsd[365];

    auto g_yy_0_0_0_y_xy_0_xx = buffer_2000_pdsd[366];

    auto g_yy_0_0_0_y_xy_0_xy = buffer_2000_pdsd[367];

    auto g_yy_0_0_0_y_xy_0_xz = buffer_2000_pdsd[368];

    auto g_yy_0_0_0_y_xy_0_yy = buffer_2000_pdsd[369];

    auto g_yy_0_0_0_y_xy_0_yz = buffer_2000_pdsd[370];

    auto g_yy_0_0_0_y_xy_0_zz = buffer_2000_pdsd[371];

    auto g_yy_0_0_0_y_xz_0_xx = buffer_2000_pdsd[372];

    auto g_yy_0_0_0_y_xz_0_xy = buffer_2000_pdsd[373];

    auto g_yy_0_0_0_y_xz_0_xz = buffer_2000_pdsd[374];

    auto g_yy_0_0_0_y_xz_0_yy = buffer_2000_pdsd[375];

    auto g_yy_0_0_0_y_xz_0_yz = buffer_2000_pdsd[376];

    auto g_yy_0_0_0_y_xz_0_zz = buffer_2000_pdsd[377];

    auto g_yy_0_0_0_y_yy_0_xx = buffer_2000_pdsd[378];

    auto g_yy_0_0_0_y_yy_0_xy = buffer_2000_pdsd[379];

    auto g_yy_0_0_0_y_yy_0_xz = buffer_2000_pdsd[380];

    auto g_yy_0_0_0_y_yy_0_yy = buffer_2000_pdsd[381];

    auto g_yy_0_0_0_y_yy_0_yz = buffer_2000_pdsd[382];

    auto g_yy_0_0_0_y_yy_0_zz = buffer_2000_pdsd[383];

    auto g_yy_0_0_0_y_yz_0_xx = buffer_2000_pdsd[384];

    auto g_yy_0_0_0_y_yz_0_xy = buffer_2000_pdsd[385];

    auto g_yy_0_0_0_y_yz_0_xz = buffer_2000_pdsd[386];

    auto g_yy_0_0_0_y_yz_0_yy = buffer_2000_pdsd[387];

    auto g_yy_0_0_0_y_yz_0_yz = buffer_2000_pdsd[388];

    auto g_yy_0_0_0_y_yz_0_zz = buffer_2000_pdsd[389];

    auto g_yy_0_0_0_y_zz_0_xx = buffer_2000_pdsd[390];

    auto g_yy_0_0_0_y_zz_0_xy = buffer_2000_pdsd[391];

    auto g_yy_0_0_0_y_zz_0_xz = buffer_2000_pdsd[392];

    auto g_yy_0_0_0_y_zz_0_yy = buffer_2000_pdsd[393];

    auto g_yy_0_0_0_y_zz_0_yz = buffer_2000_pdsd[394];

    auto g_yy_0_0_0_y_zz_0_zz = buffer_2000_pdsd[395];

    auto g_yy_0_0_0_z_xx_0_xx = buffer_2000_pdsd[396];

    auto g_yy_0_0_0_z_xx_0_xy = buffer_2000_pdsd[397];

    auto g_yy_0_0_0_z_xx_0_xz = buffer_2000_pdsd[398];

    auto g_yy_0_0_0_z_xx_0_yy = buffer_2000_pdsd[399];

    auto g_yy_0_0_0_z_xx_0_yz = buffer_2000_pdsd[400];

    auto g_yy_0_0_0_z_xx_0_zz = buffer_2000_pdsd[401];

    auto g_yy_0_0_0_z_xy_0_xx = buffer_2000_pdsd[402];

    auto g_yy_0_0_0_z_xy_0_xy = buffer_2000_pdsd[403];

    auto g_yy_0_0_0_z_xy_0_xz = buffer_2000_pdsd[404];

    auto g_yy_0_0_0_z_xy_0_yy = buffer_2000_pdsd[405];

    auto g_yy_0_0_0_z_xy_0_yz = buffer_2000_pdsd[406];

    auto g_yy_0_0_0_z_xy_0_zz = buffer_2000_pdsd[407];

    auto g_yy_0_0_0_z_xz_0_xx = buffer_2000_pdsd[408];

    auto g_yy_0_0_0_z_xz_0_xy = buffer_2000_pdsd[409];

    auto g_yy_0_0_0_z_xz_0_xz = buffer_2000_pdsd[410];

    auto g_yy_0_0_0_z_xz_0_yy = buffer_2000_pdsd[411];

    auto g_yy_0_0_0_z_xz_0_yz = buffer_2000_pdsd[412];

    auto g_yy_0_0_0_z_xz_0_zz = buffer_2000_pdsd[413];

    auto g_yy_0_0_0_z_yy_0_xx = buffer_2000_pdsd[414];

    auto g_yy_0_0_0_z_yy_0_xy = buffer_2000_pdsd[415];

    auto g_yy_0_0_0_z_yy_0_xz = buffer_2000_pdsd[416];

    auto g_yy_0_0_0_z_yy_0_yy = buffer_2000_pdsd[417];

    auto g_yy_0_0_0_z_yy_0_yz = buffer_2000_pdsd[418];

    auto g_yy_0_0_0_z_yy_0_zz = buffer_2000_pdsd[419];

    auto g_yy_0_0_0_z_yz_0_xx = buffer_2000_pdsd[420];

    auto g_yy_0_0_0_z_yz_0_xy = buffer_2000_pdsd[421];

    auto g_yy_0_0_0_z_yz_0_xz = buffer_2000_pdsd[422];

    auto g_yy_0_0_0_z_yz_0_yy = buffer_2000_pdsd[423];

    auto g_yy_0_0_0_z_yz_0_yz = buffer_2000_pdsd[424];

    auto g_yy_0_0_0_z_yz_0_zz = buffer_2000_pdsd[425];

    auto g_yy_0_0_0_z_zz_0_xx = buffer_2000_pdsd[426];

    auto g_yy_0_0_0_z_zz_0_xy = buffer_2000_pdsd[427];

    auto g_yy_0_0_0_z_zz_0_xz = buffer_2000_pdsd[428];

    auto g_yy_0_0_0_z_zz_0_yy = buffer_2000_pdsd[429];

    auto g_yy_0_0_0_z_zz_0_yz = buffer_2000_pdsd[430];

    auto g_yy_0_0_0_z_zz_0_zz = buffer_2000_pdsd[431];

    auto g_yz_0_0_0_x_xx_0_xx = buffer_2000_pdsd[432];

    auto g_yz_0_0_0_x_xx_0_xy = buffer_2000_pdsd[433];

    auto g_yz_0_0_0_x_xx_0_xz = buffer_2000_pdsd[434];

    auto g_yz_0_0_0_x_xx_0_yy = buffer_2000_pdsd[435];

    auto g_yz_0_0_0_x_xx_0_yz = buffer_2000_pdsd[436];

    auto g_yz_0_0_0_x_xx_0_zz = buffer_2000_pdsd[437];

    auto g_yz_0_0_0_x_xy_0_xx = buffer_2000_pdsd[438];

    auto g_yz_0_0_0_x_xy_0_xy = buffer_2000_pdsd[439];

    auto g_yz_0_0_0_x_xy_0_xz = buffer_2000_pdsd[440];

    auto g_yz_0_0_0_x_xy_0_yy = buffer_2000_pdsd[441];

    auto g_yz_0_0_0_x_xy_0_yz = buffer_2000_pdsd[442];

    auto g_yz_0_0_0_x_xy_0_zz = buffer_2000_pdsd[443];

    auto g_yz_0_0_0_x_xz_0_xx = buffer_2000_pdsd[444];

    auto g_yz_0_0_0_x_xz_0_xy = buffer_2000_pdsd[445];

    auto g_yz_0_0_0_x_xz_0_xz = buffer_2000_pdsd[446];

    auto g_yz_0_0_0_x_xz_0_yy = buffer_2000_pdsd[447];

    auto g_yz_0_0_0_x_xz_0_yz = buffer_2000_pdsd[448];

    auto g_yz_0_0_0_x_xz_0_zz = buffer_2000_pdsd[449];

    auto g_yz_0_0_0_x_yy_0_xx = buffer_2000_pdsd[450];

    auto g_yz_0_0_0_x_yy_0_xy = buffer_2000_pdsd[451];

    auto g_yz_0_0_0_x_yy_0_xz = buffer_2000_pdsd[452];

    auto g_yz_0_0_0_x_yy_0_yy = buffer_2000_pdsd[453];

    auto g_yz_0_0_0_x_yy_0_yz = buffer_2000_pdsd[454];

    auto g_yz_0_0_0_x_yy_0_zz = buffer_2000_pdsd[455];

    auto g_yz_0_0_0_x_yz_0_xx = buffer_2000_pdsd[456];

    auto g_yz_0_0_0_x_yz_0_xy = buffer_2000_pdsd[457];

    auto g_yz_0_0_0_x_yz_0_xz = buffer_2000_pdsd[458];

    auto g_yz_0_0_0_x_yz_0_yy = buffer_2000_pdsd[459];

    auto g_yz_0_0_0_x_yz_0_yz = buffer_2000_pdsd[460];

    auto g_yz_0_0_0_x_yz_0_zz = buffer_2000_pdsd[461];

    auto g_yz_0_0_0_x_zz_0_xx = buffer_2000_pdsd[462];

    auto g_yz_0_0_0_x_zz_0_xy = buffer_2000_pdsd[463];

    auto g_yz_0_0_0_x_zz_0_xz = buffer_2000_pdsd[464];

    auto g_yz_0_0_0_x_zz_0_yy = buffer_2000_pdsd[465];

    auto g_yz_0_0_0_x_zz_0_yz = buffer_2000_pdsd[466];

    auto g_yz_0_0_0_x_zz_0_zz = buffer_2000_pdsd[467];

    auto g_yz_0_0_0_y_xx_0_xx = buffer_2000_pdsd[468];

    auto g_yz_0_0_0_y_xx_0_xy = buffer_2000_pdsd[469];

    auto g_yz_0_0_0_y_xx_0_xz = buffer_2000_pdsd[470];

    auto g_yz_0_0_0_y_xx_0_yy = buffer_2000_pdsd[471];

    auto g_yz_0_0_0_y_xx_0_yz = buffer_2000_pdsd[472];

    auto g_yz_0_0_0_y_xx_0_zz = buffer_2000_pdsd[473];

    auto g_yz_0_0_0_y_xy_0_xx = buffer_2000_pdsd[474];

    auto g_yz_0_0_0_y_xy_0_xy = buffer_2000_pdsd[475];

    auto g_yz_0_0_0_y_xy_0_xz = buffer_2000_pdsd[476];

    auto g_yz_0_0_0_y_xy_0_yy = buffer_2000_pdsd[477];

    auto g_yz_0_0_0_y_xy_0_yz = buffer_2000_pdsd[478];

    auto g_yz_0_0_0_y_xy_0_zz = buffer_2000_pdsd[479];

    auto g_yz_0_0_0_y_xz_0_xx = buffer_2000_pdsd[480];

    auto g_yz_0_0_0_y_xz_0_xy = buffer_2000_pdsd[481];

    auto g_yz_0_0_0_y_xz_0_xz = buffer_2000_pdsd[482];

    auto g_yz_0_0_0_y_xz_0_yy = buffer_2000_pdsd[483];

    auto g_yz_0_0_0_y_xz_0_yz = buffer_2000_pdsd[484];

    auto g_yz_0_0_0_y_xz_0_zz = buffer_2000_pdsd[485];

    auto g_yz_0_0_0_y_yy_0_xx = buffer_2000_pdsd[486];

    auto g_yz_0_0_0_y_yy_0_xy = buffer_2000_pdsd[487];

    auto g_yz_0_0_0_y_yy_0_xz = buffer_2000_pdsd[488];

    auto g_yz_0_0_0_y_yy_0_yy = buffer_2000_pdsd[489];

    auto g_yz_0_0_0_y_yy_0_yz = buffer_2000_pdsd[490];

    auto g_yz_0_0_0_y_yy_0_zz = buffer_2000_pdsd[491];

    auto g_yz_0_0_0_y_yz_0_xx = buffer_2000_pdsd[492];

    auto g_yz_0_0_0_y_yz_0_xy = buffer_2000_pdsd[493];

    auto g_yz_0_0_0_y_yz_0_xz = buffer_2000_pdsd[494];

    auto g_yz_0_0_0_y_yz_0_yy = buffer_2000_pdsd[495];

    auto g_yz_0_0_0_y_yz_0_yz = buffer_2000_pdsd[496];

    auto g_yz_0_0_0_y_yz_0_zz = buffer_2000_pdsd[497];

    auto g_yz_0_0_0_y_zz_0_xx = buffer_2000_pdsd[498];

    auto g_yz_0_0_0_y_zz_0_xy = buffer_2000_pdsd[499];

    auto g_yz_0_0_0_y_zz_0_xz = buffer_2000_pdsd[500];

    auto g_yz_0_0_0_y_zz_0_yy = buffer_2000_pdsd[501];

    auto g_yz_0_0_0_y_zz_0_yz = buffer_2000_pdsd[502];

    auto g_yz_0_0_0_y_zz_0_zz = buffer_2000_pdsd[503];

    auto g_yz_0_0_0_z_xx_0_xx = buffer_2000_pdsd[504];

    auto g_yz_0_0_0_z_xx_0_xy = buffer_2000_pdsd[505];

    auto g_yz_0_0_0_z_xx_0_xz = buffer_2000_pdsd[506];

    auto g_yz_0_0_0_z_xx_0_yy = buffer_2000_pdsd[507];

    auto g_yz_0_0_0_z_xx_0_yz = buffer_2000_pdsd[508];

    auto g_yz_0_0_0_z_xx_0_zz = buffer_2000_pdsd[509];

    auto g_yz_0_0_0_z_xy_0_xx = buffer_2000_pdsd[510];

    auto g_yz_0_0_0_z_xy_0_xy = buffer_2000_pdsd[511];

    auto g_yz_0_0_0_z_xy_0_xz = buffer_2000_pdsd[512];

    auto g_yz_0_0_0_z_xy_0_yy = buffer_2000_pdsd[513];

    auto g_yz_0_0_0_z_xy_0_yz = buffer_2000_pdsd[514];

    auto g_yz_0_0_0_z_xy_0_zz = buffer_2000_pdsd[515];

    auto g_yz_0_0_0_z_xz_0_xx = buffer_2000_pdsd[516];

    auto g_yz_0_0_0_z_xz_0_xy = buffer_2000_pdsd[517];

    auto g_yz_0_0_0_z_xz_0_xz = buffer_2000_pdsd[518];

    auto g_yz_0_0_0_z_xz_0_yy = buffer_2000_pdsd[519];

    auto g_yz_0_0_0_z_xz_0_yz = buffer_2000_pdsd[520];

    auto g_yz_0_0_0_z_xz_0_zz = buffer_2000_pdsd[521];

    auto g_yz_0_0_0_z_yy_0_xx = buffer_2000_pdsd[522];

    auto g_yz_0_0_0_z_yy_0_xy = buffer_2000_pdsd[523];

    auto g_yz_0_0_0_z_yy_0_xz = buffer_2000_pdsd[524];

    auto g_yz_0_0_0_z_yy_0_yy = buffer_2000_pdsd[525];

    auto g_yz_0_0_0_z_yy_0_yz = buffer_2000_pdsd[526];

    auto g_yz_0_0_0_z_yy_0_zz = buffer_2000_pdsd[527];

    auto g_yz_0_0_0_z_yz_0_xx = buffer_2000_pdsd[528];

    auto g_yz_0_0_0_z_yz_0_xy = buffer_2000_pdsd[529];

    auto g_yz_0_0_0_z_yz_0_xz = buffer_2000_pdsd[530];

    auto g_yz_0_0_0_z_yz_0_yy = buffer_2000_pdsd[531];

    auto g_yz_0_0_0_z_yz_0_yz = buffer_2000_pdsd[532];

    auto g_yz_0_0_0_z_yz_0_zz = buffer_2000_pdsd[533];

    auto g_yz_0_0_0_z_zz_0_xx = buffer_2000_pdsd[534];

    auto g_yz_0_0_0_z_zz_0_xy = buffer_2000_pdsd[535];

    auto g_yz_0_0_0_z_zz_0_xz = buffer_2000_pdsd[536];

    auto g_yz_0_0_0_z_zz_0_yy = buffer_2000_pdsd[537];

    auto g_yz_0_0_0_z_zz_0_yz = buffer_2000_pdsd[538];

    auto g_yz_0_0_0_z_zz_0_zz = buffer_2000_pdsd[539];

    auto g_zz_0_0_0_x_xx_0_xx = buffer_2000_pdsd[540];

    auto g_zz_0_0_0_x_xx_0_xy = buffer_2000_pdsd[541];

    auto g_zz_0_0_0_x_xx_0_xz = buffer_2000_pdsd[542];

    auto g_zz_0_0_0_x_xx_0_yy = buffer_2000_pdsd[543];

    auto g_zz_0_0_0_x_xx_0_yz = buffer_2000_pdsd[544];

    auto g_zz_0_0_0_x_xx_0_zz = buffer_2000_pdsd[545];

    auto g_zz_0_0_0_x_xy_0_xx = buffer_2000_pdsd[546];

    auto g_zz_0_0_0_x_xy_0_xy = buffer_2000_pdsd[547];

    auto g_zz_0_0_0_x_xy_0_xz = buffer_2000_pdsd[548];

    auto g_zz_0_0_0_x_xy_0_yy = buffer_2000_pdsd[549];

    auto g_zz_0_0_0_x_xy_0_yz = buffer_2000_pdsd[550];

    auto g_zz_0_0_0_x_xy_0_zz = buffer_2000_pdsd[551];

    auto g_zz_0_0_0_x_xz_0_xx = buffer_2000_pdsd[552];

    auto g_zz_0_0_0_x_xz_0_xy = buffer_2000_pdsd[553];

    auto g_zz_0_0_0_x_xz_0_xz = buffer_2000_pdsd[554];

    auto g_zz_0_0_0_x_xz_0_yy = buffer_2000_pdsd[555];

    auto g_zz_0_0_0_x_xz_0_yz = buffer_2000_pdsd[556];

    auto g_zz_0_0_0_x_xz_0_zz = buffer_2000_pdsd[557];

    auto g_zz_0_0_0_x_yy_0_xx = buffer_2000_pdsd[558];

    auto g_zz_0_0_0_x_yy_0_xy = buffer_2000_pdsd[559];

    auto g_zz_0_0_0_x_yy_0_xz = buffer_2000_pdsd[560];

    auto g_zz_0_0_0_x_yy_0_yy = buffer_2000_pdsd[561];

    auto g_zz_0_0_0_x_yy_0_yz = buffer_2000_pdsd[562];

    auto g_zz_0_0_0_x_yy_0_zz = buffer_2000_pdsd[563];

    auto g_zz_0_0_0_x_yz_0_xx = buffer_2000_pdsd[564];

    auto g_zz_0_0_0_x_yz_0_xy = buffer_2000_pdsd[565];

    auto g_zz_0_0_0_x_yz_0_xz = buffer_2000_pdsd[566];

    auto g_zz_0_0_0_x_yz_0_yy = buffer_2000_pdsd[567];

    auto g_zz_0_0_0_x_yz_0_yz = buffer_2000_pdsd[568];

    auto g_zz_0_0_0_x_yz_0_zz = buffer_2000_pdsd[569];

    auto g_zz_0_0_0_x_zz_0_xx = buffer_2000_pdsd[570];

    auto g_zz_0_0_0_x_zz_0_xy = buffer_2000_pdsd[571];

    auto g_zz_0_0_0_x_zz_0_xz = buffer_2000_pdsd[572];

    auto g_zz_0_0_0_x_zz_0_yy = buffer_2000_pdsd[573];

    auto g_zz_0_0_0_x_zz_0_yz = buffer_2000_pdsd[574];

    auto g_zz_0_0_0_x_zz_0_zz = buffer_2000_pdsd[575];

    auto g_zz_0_0_0_y_xx_0_xx = buffer_2000_pdsd[576];

    auto g_zz_0_0_0_y_xx_0_xy = buffer_2000_pdsd[577];

    auto g_zz_0_0_0_y_xx_0_xz = buffer_2000_pdsd[578];

    auto g_zz_0_0_0_y_xx_0_yy = buffer_2000_pdsd[579];

    auto g_zz_0_0_0_y_xx_0_yz = buffer_2000_pdsd[580];

    auto g_zz_0_0_0_y_xx_0_zz = buffer_2000_pdsd[581];

    auto g_zz_0_0_0_y_xy_0_xx = buffer_2000_pdsd[582];

    auto g_zz_0_0_0_y_xy_0_xy = buffer_2000_pdsd[583];

    auto g_zz_0_0_0_y_xy_0_xz = buffer_2000_pdsd[584];

    auto g_zz_0_0_0_y_xy_0_yy = buffer_2000_pdsd[585];

    auto g_zz_0_0_0_y_xy_0_yz = buffer_2000_pdsd[586];

    auto g_zz_0_0_0_y_xy_0_zz = buffer_2000_pdsd[587];

    auto g_zz_0_0_0_y_xz_0_xx = buffer_2000_pdsd[588];

    auto g_zz_0_0_0_y_xz_0_xy = buffer_2000_pdsd[589];

    auto g_zz_0_0_0_y_xz_0_xz = buffer_2000_pdsd[590];

    auto g_zz_0_0_0_y_xz_0_yy = buffer_2000_pdsd[591];

    auto g_zz_0_0_0_y_xz_0_yz = buffer_2000_pdsd[592];

    auto g_zz_0_0_0_y_xz_0_zz = buffer_2000_pdsd[593];

    auto g_zz_0_0_0_y_yy_0_xx = buffer_2000_pdsd[594];

    auto g_zz_0_0_0_y_yy_0_xy = buffer_2000_pdsd[595];

    auto g_zz_0_0_0_y_yy_0_xz = buffer_2000_pdsd[596];

    auto g_zz_0_0_0_y_yy_0_yy = buffer_2000_pdsd[597];

    auto g_zz_0_0_0_y_yy_0_yz = buffer_2000_pdsd[598];

    auto g_zz_0_0_0_y_yy_0_zz = buffer_2000_pdsd[599];

    auto g_zz_0_0_0_y_yz_0_xx = buffer_2000_pdsd[600];

    auto g_zz_0_0_0_y_yz_0_xy = buffer_2000_pdsd[601];

    auto g_zz_0_0_0_y_yz_0_xz = buffer_2000_pdsd[602];

    auto g_zz_0_0_0_y_yz_0_yy = buffer_2000_pdsd[603];

    auto g_zz_0_0_0_y_yz_0_yz = buffer_2000_pdsd[604];

    auto g_zz_0_0_0_y_yz_0_zz = buffer_2000_pdsd[605];

    auto g_zz_0_0_0_y_zz_0_xx = buffer_2000_pdsd[606];

    auto g_zz_0_0_0_y_zz_0_xy = buffer_2000_pdsd[607];

    auto g_zz_0_0_0_y_zz_0_xz = buffer_2000_pdsd[608];

    auto g_zz_0_0_0_y_zz_0_yy = buffer_2000_pdsd[609];

    auto g_zz_0_0_0_y_zz_0_yz = buffer_2000_pdsd[610];

    auto g_zz_0_0_0_y_zz_0_zz = buffer_2000_pdsd[611];

    auto g_zz_0_0_0_z_xx_0_xx = buffer_2000_pdsd[612];

    auto g_zz_0_0_0_z_xx_0_xy = buffer_2000_pdsd[613];

    auto g_zz_0_0_0_z_xx_0_xz = buffer_2000_pdsd[614];

    auto g_zz_0_0_0_z_xx_0_yy = buffer_2000_pdsd[615];

    auto g_zz_0_0_0_z_xx_0_yz = buffer_2000_pdsd[616];

    auto g_zz_0_0_0_z_xx_0_zz = buffer_2000_pdsd[617];

    auto g_zz_0_0_0_z_xy_0_xx = buffer_2000_pdsd[618];

    auto g_zz_0_0_0_z_xy_0_xy = buffer_2000_pdsd[619];

    auto g_zz_0_0_0_z_xy_0_xz = buffer_2000_pdsd[620];

    auto g_zz_0_0_0_z_xy_0_yy = buffer_2000_pdsd[621];

    auto g_zz_0_0_0_z_xy_0_yz = buffer_2000_pdsd[622];

    auto g_zz_0_0_0_z_xy_0_zz = buffer_2000_pdsd[623];

    auto g_zz_0_0_0_z_xz_0_xx = buffer_2000_pdsd[624];

    auto g_zz_0_0_0_z_xz_0_xy = buffer_2000_pdsd[625];

    auto g_zz_0_0_0_z_xz_0_xz = buffer_2000_pdsd[626];

    auto g_zz_0_0_0_z_xz_0_yy = buffer_2000_pdsd[627];

    auto g_zz_0_0_0_z_xz_0_yz = buffer_2000_pdsd[628];

    auto g_zz_0_0_0_z_xz_0_zz = buffer_2000_pdsd[629];

    auto g_zz_0_0_0_z_yy_0_xx = buffer_2000_pdsd[630];

    auto g_zz_0_0_0_z_yy_0_xy = buffer_2000_pdsd[631];

    auto g_zz_0_0_0_z_yy_0_xz = buffer_2000_pdsd[632];

    auto g_zz_0_0_0_z_yy_0_yy = buffer_2000_pdsd[633];

    auto g_zz_0_0_0_z_yy_0_yz = buffer_2000_pdsd[634];

    auto g_zz_0_0_0_z_yy_0_zz = buffer_2000_pdsd[635];

    auto g_zz_0_0_0_z_yz_0_xx = buffer_2000_pdsd[636];

    auto g_zz_0_0_0_z_yz_0_xy = buffer_2000_pdsd[637];

    auto g_zz_0_0_0_z_yz_0_xz = buffer_2000_pdsd[638];

    auto g_zz_0_0_0_z_yz_0_yy = buffer_2000_pdsd[639];

    auto g_zz_0_0_0_z_yz_0_yz = buffer_2000_pdsd[640];

    auto g_zz_0_0_0_z_yz_0_zz = buffer_2000_pdsd[641];

    auto g_zz_0_0_0_z_zz_0_xx = buffer_2000_pdsd[642];

    auto g_zz_0_0_0_z_zz_0_xy = buffer_2000_pdsd[643];

    auto g_zz_0_0_0_z_zz_0_xz = buffer_2000_pdsd[644];

    auto g_zz_0_0_0_z_zz_0_yy = buffer_2000_pdsd[645];

    auto g_zz_0_0_0_z_zz_0_yz = buffer_2000_pdsd[646];

    auto g_zz_0_0_0_z_zz_0_zz = buffer_2000_pdsd[647];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_xx_0_0_0_x_xx_0_xx, g_xx_0_0_0_x_xx_0_xy, g_xx_0_0_0_x_xx_0_xz, g_xx_0_0_0_x_xx_0_yy, g_xx_0_0_0_x_xx_0_yz, g_xx_0_0_0_x_xx_0_zz, g_xxx_xx_0_xx, g_xxx_xx_0_xy, g_xxx_xx_0_xz, g_xxx_xx_0_yy, g_xxx_xx_0_yz, g_xxx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xx_0_xx[i] = -6.0 * g_x_xx_0_xx[i] * a_exp + 4.0 * g_xxx_xx_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_0_xy[i] = -6.0 * g_x_xx_0_xy[i] * a_exp + 4.0 * g_xxx_xx_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_0_xz[i] = -6.0 * g_x_xx_0_xz[i] * a_exp + 4.0 * g_xxx_xx_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_0_yy[i] = -6.0 * g_x_xx_0_yy[i] * a_exp + 4.0 * g_xxx_xx_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_0_yz[i] = -6.0 * g_x_xx_0_yz[i] * a_exp + 4.0 * g_xxx_xx_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_0_zz[i] = -6.0 * g_x_xx_0_zz[i] * a_exp + 4.0 * g_xxx_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_xx_0_0_0_x_xy_0_xx, g_xx_0_0_0_x_xy_0_xy, g_xx_0_0_0_x_xy_0_xz, g_xx_0_0_0_x_xy_0_yy, g_xx_0_0_0_x_xy_0_yz, g_xx_0_0_0_x_xy_0_zz, g_xxx_xy_0_xx, g_xxx_xy_0_xy, g_xxx_xy_0_xz, g_xxx_xy_0_yy, g_xxx_xy_0_yz, g_xxx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xy_0_xx[i] = -6.0 * g_x_xy_0_xx[i] * a_exp + 4.0 * g_xxx_xy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_0_xy[i] = -6.0 * g_x_xy_0_xy[i] * a_exp + 4.0 * g_xxx_xy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_0_xz[i] = -6.0 * g_x_xy_0_xz[i] * a_exp + 4.0 * g_xxx_xy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_0_yy[i] = -6.0 * g_x_xy_0_yy[i] * a_exp + 4.0 * g_xxx_xy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_0_yz[i] = -6.0 * g_x_xy_0_yz[i] * a_exp + 4.0 * g_xxx_xy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_0_zz[i] = -6.0 * g_x_xy_0_zz[i] * a_exp + 4.0 * g_xxx_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_xx_0_0_0_x_xz_0_xx, g_xx_0_0_0_x_xz_0_xy, g_xx_0_0_0_x_xz_0_xz, g_xx_0_0_0_x_xz_0_yy, g_xx_0_0_0_x_xz_0_yz, g_xx_0_0_0_x_xz_0_zz, g_xxx_xz_0_xx, g_xxx_xz_0_xy, g_xxx_xz_0_xz, g_xxx_xz_0_yy, g_xxx_xz_0_yz, g_xxx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xz_0_xx[i] = -6.0 * g_x_xz_0_xx[i] * a_exp + 4.0 * g_xxx_xz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_0_xy[i] = -6.0 * g_x_xz_0_xy[i] * a_exp + 4.0 * g_xxx_xz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_0_xz[i] = -6.0 * g_x_xz_0_xz[i] * a_exp + 4.0 * g_xxx_xz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_0_yy[i] = -6.0 * g_x_xz_0_yy[i] * a_exp + 4.0 * g_xxx_xz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_0_yz[i] = -6.0 * g_x_xz_0_yz[i] * a_exp + 4.0 * g_xxx_xz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_0_zz[i] = -6.0 * g_x_xz_0_zz[i] * a_exp + 4.0 * g_xxx_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_xx_0_0_0_x_yy_0_xx, g_xx_0_0_0_x_yy_0_xy, g_xx_0_0_0_x_yy_0_xz, g_xx_0_0_0_x_yy_0_yy, g_xx_0_0_0_x_yy_0_yz, g_xx_0_0_0_x_yy_0_zz, g_xxx_yy_0_xx, g_xxx_yy_0_xy, g_xxx_yy_0_xz, g_xxx_yy_0_yy, g_xxx_yy_0_yz, g_xxx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_yy_0_xx[i] = -6.0 * g_x_yy_0_xx[i] * a_exp + 4.0 * g_xxx_yy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_0_xy[i] = -6.0 * g_x_yy_0_xy[i] * a_exp + 4.0 * g_xxx_yy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_0_xz[i] = -6.0 * g_x_yy_0_xz[i] * a_exp + 4.0 * g_xxx_yy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_0_yy[i] = -6.0 * g_x_yy_0_yy[i] * a_exp + 4.0 * g_xxx_yy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_0_yz[i] = -6.0 * g_x_yy_0_yz[i] * a_exp + 4.0 * g_xxx_yy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_0_zz[i] = -6.0 * g_x_yy_0_zz[i] * a_exp + 4.0 * g_xxx_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_xx_0_0_0_x_yz_0_xx, g_xx_0_0_0_x_yz_0_xy, g_xx_0_0_0_x_yz_0_xz, g_xx_0_0_0_x_yz_0_yy, g_xx_0_0_0_x_yz_0_yz, g_xx_0_0_0_x_yz_0_zz, g_xxx_yz_0_xx, g_xxx_yz_0_xy, g_xxx_yz_0_xz, g_xxx_yz_0_yy, g_xxx_yz_0_yz, g_xxx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_yz_0_xx[i] = -6.0 * g_x_yz_0_xx[i] * a_exp + 4.0 * g_xxx_yz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_0_xy[i] = -6.0 * g_x_yz_0_xy[i] * a_exp + 4.0 * g_xxx_yz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_0_xz[i] = -6.0 * g_x_yz_0_xz[i] * a_exp + 4.0 * g_xxx_yz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_0_yy[i] = -6.0 * g_x_yz_0_yy[i] * a_exp + 4.0 * g_xxx_yz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_0_yz[i] = -6.0 * g_x_yz_0_yz[i] * a_exp + 4.0 * g_xxx_yz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_0_zz[i] = -6.0 * g_x_yz_0_zz[i] * a_exp + 4.0 * g_xxx_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_xx_0_0_0_x_zz_0_xx, g_xx_0_0_0_x_zz_0_xy, g_xx_0_0_0_x_zz_0_xz, g_xx_0_0_0_x_zz_0_yy, g_xx_0_0_0_x_zz_0_yz, g_xx_0_0_0_x_zz_0_zz, g_xxx_zz_0_xx, g_xxx_zz_0_xy, g_xxx_zz_0_xz, g_xxx_zz_0_yy, g_xxx_zz_0_yz, g_xxx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_zz_0_xx[i] = -6.0 * g_x_zz_0_xx[i] * a_exp + 4.0 * g_xxx_zz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_0_xy[i] = -6.0 * g_x_zz_0_xy[i] * a_exp + 4.0 * g_xxx_zz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_0_xz[i] = -6.0 * g_x_zz_0_xz[i] * a_exp + 4.0 * g_xxx_zz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_0_yy[i] = -6.0 * g_x_zz_0_yy[i] * a_exp + 4.0 * g_xxx_zz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_0_yz[i] = -6.0 * g_x_zz_0_yz[i] * a_exp + 4.0 * g_xxx_zz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_0_zz[i] = -6.0 * g_x_zz_0_zz[i] * a_exp + 4.0 * g_xxx_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xx_0_0_0_y_xx_0_xx, g_xx_0_0_0_y_xx_0_xy, g_xx_0_0_0_y_xx_0_xz, g_xx_0_0_0_y_xx_0_yy, g_xx_0_0_0_y_xx_0_yz, g_xx_0_0_0_y_xx_0_zz, g_xxy_xx_0_xx, g_xxy_xx_0_xy, g_xxy_xx_0_xz, g_xxy_xx_0_yy, g_xxy_xx_0_yz, g_xxy_xx_0_zz, g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xx_0_xx[i] = -2.0 * g_y_xx_0_xx[i] * a_exp + 4.0 * g_xxy_xx_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_0_xy[i] = -2.0 * g_y_xx_0_xy[i] * a_exp + 4.0 * g_xxy_xx_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_0_xz[i] = -2.0 * g_y_xx_0_xz[i] * a_exp + 4.0 * g_xxy_xx_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_0_yy[i] = -2.0 * g_y_xx_0_yy[i] * a_exp + 4.0 * g_xxy_xx_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_0_yz[i] = -2.0 * g_y_xx_0_yz[i] * a_exp + 4.0 * g_xxy_xx_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_0_zz[i] = -2.0 * g_y_xx_0_zz[i] * a_exp + 4.0 * g_xxy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_xx_0_0_0_y_xy_0_xx, g_xx_0_0_0_y_xy_0_xy, g_xx_0_0_0_y_xy_0_xz, g_xx_0_0_0_y_xy_0_yy, g_xx_0_0_0_y_xy_0_yz, g_xx_0_0_0_y_xy_0_zz, g_xxy_xy_0_xx, g_xxy_xy_0_xy, g_xxy_xy_0_xz, g_xxy_xy_0_yy, g_xxy_xy_0_yz, g_xxy_xy_0_zz, g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xy_0_xx[i] = -2.0 * g_y_xy_0_xx[i] * a_exp + 4.0 * g_xxy_xy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_0_xy[i] = -2.0 * g_y_xy_0_xy[i] * a_exp + 4.0 * g_xxy_xy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_0_xz[i] = -2.0 * g_y_xy_0_xz[i] * a_exp + 4.0 * g_xxy_xy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_0_yy[i] = -2.0 * g_y_xy_0_yy[i] * a_exp + 4.0 * g_xxy_xy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_0_yz[i] = -2.0 * g_y_xy_0_yz[i] * a_exp + 4.0 * g_xxy_xy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_0_zz[i] = -2.0 * g_y_xy_0_zz[i] * a_exp + 4.0 * g_xxy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_xx_0_0_0_y_xz_0_xx, g_xx_0_0_0_y_xz_0_xy, g_xx_0_0_0_y_xz_0_xz, g_xx_0_0_0_y_xz_0_yy, g_xx_0_0_0_y_xz_0_yz, g_xx_0_0_0_y_xz_0_zz, g_xxy_xz_0_xx, g_xxy_xz_0_xy, g_xxy_xz_0_xz, g_xxy_xz_0_yy, g_xxy_xz_0_yz, g_xxy_xz_0_zz, g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xz_0_xx[i] = -2.0 * g_y_xz_0_xx[i] * a_exp + 4.0 * g_xxy_xz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_0_xy[i] = -2.0 * g_y_xz_0_xy[i] * a_exp + 4.0 * g_xxy_xz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_0_xz[i] = -2.0 * g_y_xz_0_xz[i] * a_exp + 4.0 * g_xxy_xz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_0_yy[i] = -2.0 * g_y_xz_0_yy[i] * a_exp + 4.0 * g_xxy_xz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_0_yz[i] = -2.0 * g_y_xz_0_yz[i] * a_exp + 4.0 * g_xxy_xz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_0_zz[i] = -2.0 * g_y_xz_0_zz[i] * a_exp + 4.0 * g_xxy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xx_0_0_0_y_yy_0_xx, g_xx_0_0_0_y_yy_0_xy, g_xx_0_0_0_y_yy_0_xz, g_xx_0_0_0_y_yy_0_yy, g_xx_0_0_0_y_yy_0_yz, g_xx_0_0_0_y_yy_0_zz, g_xxy_yy_0_xx, g_xxy_yy_0_xy, g_xxy_yy_0_xz, g_xxy_yy_0_yy, g_xxy_yy_0_yz, g_xxy_yy_0_zz, g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_yy_0_xx[i] = -2.0 * g_y_yy_0_xx[i] * a_exp + 4.0 * g_xxy_yy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_0_xy[i] = -2.0 * g_y_yy_0_xy[i] * a_exp + 4.0 * g_xxy_yy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_0_xz[i] = -2.0 * g_y_yy_0_xz[i] * a_exp + 4.0 * g_xxy_yy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_0_yy[i] = -2.0 * g_y_yy_0_yy[i] * a_exp + 4.0 * g_xxy_yy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_0_yz[i] = -2.0 * g_y_yy_0_yz[i] * a_exp + 4.0 * g_xxy_yy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_0_zz[i] = -2.0 * g_y_yy_0_zz[i] * a_exp + 4.0 * g_xxy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_xx_0_0_0_y_yz_0_xx, g_xx_0_0_0_y_yz_0_xy, g_xx_0_0_0_y_yz_0_xz, g_xx_0_0_0_y_yz_0_yy, g_xx_0_0_0_y_yz_0_yz, g_xx_0_0_0_y_yz_0_zz, g_xxy_yz_0_xx, g_xxy_yz_0_xy, g_xxy_yz_0_xz, g_xxy_yz_0_yy, g_xxy_yz_0_yz, g_xxy_yz_0_zz, g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_yz_0_xx[i] = -2.0 * g_y_yz_0_xx[i] * a_exp + 4.0 * g_xxy_yz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_0_xy[i] = -2.0 * g_y_yz_0_xy[i] * a_exp + 4.0 * g_xxy_yz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_0_xz[i] = -2.0 * g_y_yz_0_xz[i] * a_exp + 4.0 * g_xxy_yz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_0_yy[i] = -2.0 * g_y_yz_0_yy[i] * a_exp + 4.0 * g_xxy_yz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_0_yz[i] = -2.0 * g_y_yz_0_yz[i] * a_exp + 4.0 * g_xxy_yz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_0_zz[i] = -2.0 * g_y_yz_0_zz[i] * a_exp + 4.0 * g_xxy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_xx_0_0_0_y_zz_0_xx, g_xx_0_0_0_y_zz_0_xy, g_xx_0_0_0_y_zz_0_xz, g_xx_0_0_0_y_zz_0_yy, g_xx_0_0_0_y_zz_0_yz, g_xx_0_0_0_y_zz_0_zz, g_xxy_zz_0_xx, g_xxy_zz_0_xy, g_xxy_zz_0_xz, g_xxy_zz_0_yy, g_xxy_zz_0_yz, g_xxy_zz_0_zz, g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_zz_0_xx[i] = -2.0 * g_y_zz_0_xx[i] * a_exp + 4.0 * g_xxy_zz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_0_xy[i] = -2.0 * g_y_zz_0_xy[i] * a_exp + 4.0 * g_xxy_zz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_0_xz[i] = -2.0 * g_y_zz_0_xz[i] * a_exp + 4.0 * g_xxy_zz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_0_yy[i] = -2.0 * g_y_zz_0_yy[i] * a_exp + 4.0 * g_xxy_zz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_0_yz[i] = -2.0 * g_y_zz_0_yz[i] * a_exp + 4.0 * g_xxy_zz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_0_zz[i] = -2.0 * g_y_zz_0_zz[i] * a_exp + 4.0 * g_xxy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xx_0_0_0_z_xx_0_xx, g_xx_0_0_0_z_xx_0_xy, g_xx_0_0_0_z_xx_0_xz, g_xx_0_0_0_z_xx_0_yy, g_xx_0_0_0_z_xx_0_yz, g_xx_0_0_0_z_xx_0_zz, g_xxz_xx_0_xx, g_xxz_xx_0_xy, g_xxz_xx_0_xz, g_xxz_xx_0_yy, g_xxz_xx_0_yz, g_xxz_xx_0_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xx_0_xx[i] = -2.0 * g_z_xx_0_xx[i] * a_exp + 4.0 * g_xxz_xx_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_0_xy[i] = -2.0 * g_z_xx_0_xy[i] * a_exp + 4.0 * g_xxz_xx_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_0_xz[i] = -2.0 * g_z_xx_0_xz[i] * a_exp + 4.0 * g_xxz_xx_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_0_yy[i] = -2.0 * g_z_xx_0_yy[i] * a_exp + 4.0 * g_xxz_xx_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_0_yz[i] = -2.0 * g_z_xx_0_yz[i] * a_exp + 4.0 * g_xxz_xx_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_0_zz[i] = -2.0 * g_z_xx_0_zz[i] * a_exp + 4.0 * g_xxz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_xx_0_0_0_z_xy_0_xx, g_xx_0_0_0_z_xy_0_xy, g_xx_0_0_0_z_xy_0_xz, g_xx_0_0_0_z_xy_0_yy, g_xx_0_0_0_z_xy_0_yz, g_xx_0_0_0_z_xy_0_zz, g_xxz_xy_0_xx, g_xxz_xy_0_xy, g_xxz_xy_0_xz, g_xxz_xy_0_yy, g_xxz_xy_0_yz, g_xxz_xy_0_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xy_0_xx[i] = -2.0 * g_z_xy_0_xx[i] * a_exp + 4.0 * g_xxz_xy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_0_xy[i] = -2.0 * g_z_xy_0_xy[i] * a_exp + 4.0 * g_xxz_xy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_0_xz[i] = -2.0 * g_z_xy_0_xz[i] * a_exp + 4.0 * g_xxz_xy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_0_yy[i] = -2.0 * g_z_xy_0_yy[i] * a_exp + 4.0 * g_xxz_xy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_0_yz[i] = -2.0 * g_z_xy_0_yz[i] * a_exp + 4.0 * g_xxz_xy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_0_zz[i] = -2.0 * g_z_xy_0_zz[i] * a_exp + 4.0 * g_xxz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_xx_0_0_0_z_xz_0_xx, g_xx_0_0_0_z_xz_0_xy, g_xx_0_0_0_z_xz_0_xz, g_xx_0_0_0_z_xz_0_yy, g_xx_0_0_0_z_xz_0_yz, g_xx_0_0_0_z_xz_0_zz, g_xxz_xz_0_xx, g_xxz_xz_0_xy, g_xxz_xz_0_xz, g_xxz_xz_0_yy, g_xxz_xz_0_yz, g_xxz_xz_0_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xz_0_xx[i] = -2.0 * g_z_xz_0_xx[i] * a_exp + 4.0 * g_xxz_xz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_0_xy[i] = -2.0 * g_z_xz_0_xy[i] * a_exp + 4.0 * g_xxz_xz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_0_xz[i] = -2.0 * g_z_xz_0_xz[i] * a_exp + 4.0 * g_xxz_xz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_0_yy[i] = -2.0 * g_z_xz_0_yy[i] * a_exp + 4.0 * g_xxz_xz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_0_yz[i] = -2.0 * g_z_xz_0_yz[i] * a_exp + 4.0 * g_xxz_xz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_0_zz[i] = -2.0 * g_z_xz_0_zz[i] * a_exp + 4.0 * g_xxz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_xx_0_0_0_z_yy_0_xx, g_xx_0_0_0_z_yy_0_xy, g_xx_0_0_0_z_yy_0_xz, g_xx_0_0_0_z_yy_0_yy, g_xx_0_0_0_z_yy_0_yz, g_xx_0_0_0_z_yy_0_zz, g_xxz_yy_0_xx, g_xxz_yy_0_xy, g_xxz_yy_0_xz, g_xxz_yy_0_yy, g_xxz_yy_0_yz, g_xxz_yy_0_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_yy_0_xx[i] = -2.0 * g_z_yy_0_xx[i] * a_exp + 4.0 * g_xxz_yy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_0_xy[i] = -2.0 * g_z_yy_0_xy[i] * a_exp + 4.0 * g_xxz_yy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_0_xz[i] = -2.0 * g_z_yy_0_xz[i] * a_exp + 4.0 * g_xxz_yy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_0_yy[i] = -2.0 * g_z_yy_0_yy[i] * a_exp + 4.0 * g_xxz_yy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_0_yz[i] = -2.0 * g_z_yy_0_yz[i] * a_exp + 4.0 * g_xxz_yy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_0_zz[i] = -2.0 * g_z_yy_0_zz[i] * a_exp + 4.0 * g_xxz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_xx_0_0_0_z_yz_0_xx, g_xx_0_0_0_z_yz_0_xy, g_xx_0_0_0_z_yz_0_xz, g_xx_0_0_0_z_yz_0_yy, g_xx_0_0_0_z_yz_0_yz, g_xx_0_0_0_z_yz_0_zz, g_xxz_yz_0_xx, g_xxz_yz_0_xy, g_xxz_yz_0_xz, g_xxz_yz_0_yy, g_xxz_yz_0_yz, g_xxz_yz_0_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_yz_0_xx[i] = -2.0 * g_z_yz_0_xx[i] * a_exp + 4.0 * g_xxz_yz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_0_xy[i] = -2.0 * g_z_yz_0_xy[i] * a_exp + 4.0 * g_xxz_yz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_0_xz[i] = -2.0 * g_z_yz_0_xz[i] * a_exp + 4.0 * g_xxz_yz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_0_yy[i] = -2.0 * g_z_yz_0_yy[i] * a_exp + 4.0 * g_xxz_yz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_0_yz[i] = -2.0 * g_z_yz_0_yz[i] * a_exp + 4.0 * g_xxz_yz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_0_zz[i] = -2.0 * g_z_yz_0_zz[i] * a_exp + 4.0 * g_xxz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_xx_0_0_0_z_zz_0_xx, g_xx_0_0_0_z_zz_0_xy, g_xx_0_0_0_z_zz_0_xz, g_xx_0_0_0_z_zz_0_yy, g_xx_0_0_0_z_zz_0_yz, g_xx_0_0_0_z_zz_0_zz, g_xxz_zz_0_xx, g_xxz_zz_0_xy, g_xxz_zz_0_xz, g_xxz_zz_0_yy, g_xxz_zz_0_yz, g_xxz_zz_0_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_zz_0_xx[i] = -2.0 * g_z_zz_0_xx[i] * a_exp + 4.0 * g_xxz_zz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_0_xy[i] = -2.0 * g_z_zz_0_xy[i] * a_exp + 4.0 * g_xxz_zz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_0_xz[i] = -2.0 * g_z_zz_0_xz[i] * a_exp + 4.0 * g_xxz_zz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_0_yy[i] = -2.0 * g_z_zz_0_yy[i] * a_exp + 4.0 * g_xxz_zz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_0_yz[i] = -2.0 * g_z_zz_0_yz[i] * a_exp + 4.0 * g_xxz_zz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_0_zz[i] = -2.0 * g_z_zz_0_zz[i] * a_exp + 4.0 * g_xxz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xxy_xx_0_xx, g_xxy_xx_0_xy, g_xxy_xx_0_xz, g_xxy_xx_0_yy, g_xxy_xx_0_yz, g_xxy_xx_0_zz, g_xy_0_0_0_x_xx_0_xx, g_xy_0_0_0_x_xx_0_xy, g_xy_0_0_0_x_xx_0_xz, g_xy_0_0_0_x_xx_0_yy, g_xy_0_0_0_x_xx_0_yz, g_xy_0_0_0_x_xx_0_zz, g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xx_0_xx[i] = -2.0 * g_y_xx_0_xx[i] * a_exp + 4.0 * g_xxy_xx_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_0_xy[i] = -2.0 * g_y_xx_0_xy[i] * a_exp + 4.0 * g_xxy_xx_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_0_xz[i] = -2.0 * g_y_xx_0_xz[i] * a_exp + 4.0 * g_xxy_xx_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_0_yy[i] = -2.0 * g_y_xx_0_yy[i] * a_exp + 4.0 * g_xxy_xx_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_0_yz[i] = -2.0 * g_y_xx_0_yz[i] * a_exp + 4.0 * g_xxy_xx_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_0_zz[i] = -2.0 * g_y_xx_0_zz[i] * a_exp + 4.0 * g_xxy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xxy_xy_0_xx, g_xxy_xy_0_xy, g_xxy_xy_0_xz, g_xxy_xy_0_yy, g_xxy_xy_0_yz, g_xxy_xy_0_zz, g_xy_0_0_0_x_xy_0_xx, g_xy_0_0_0_x_xy_0_xy, g_xy_0_0_0_x_xy_0_xz, g_xy_0_0_0_x_xy_0_yy, g_xy_0_0_0_x_xy_0_yz, g_xy_0_0_0_x_xy_0_zz, g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xy_0_xx[i] = -2.0 * g_y_xy_0_xx[i] * a_exp + 4.0 * g_xxy_xy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_0_xy[i] = -2.0 * g_y_xy_0_xy[i] * a_exp + 4.0 * g_xxy_xy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_0_xz[i] = -2.0 * g_y_xy_0_xz[i] * a_exp + 4.0 * g_xxy_xy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_0_yy[i] = -2.0 * g_y_xy_0_yy[i] * a_exp + 4.0 * g_xxy_xy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_0_yz[i] = -2.0 * g_y_xy_0_yz[i] * a_exp + 4.0 * g_xxy_xy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_0_zz[i] = -2.0 * g_y_xy_0_zz[i] * a_exp + 4.0 * g_xxy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xxy_xz_0_xx, g_xxy_xz_0_xy, g_xxy_xz_0_xz, g_xxy_xz_0_yy, g_xxy_xz_0_yz, g_xxy_xz_0_zz, g_xy_0_0_0_x_xz_0_xx, g_xy_0_0_0_x_xz_0_xy, g_xy_0_0_0_x_xz_0_xz, g_xy_0_0_0_x_xz_0_yy, g_xy_0_0_0_x_xz_0_yz, g_xy_0_0_0_x_xz_0_zz, g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xz_0_xx[i] = -2.0 * g_y_xz_0_xx[i] * a_exp + 4.0 * g_xxy_xz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_0_xy[i] = -2.0 * g_y_xz_0_xy[i] * a_exp + 4.0 * g_xxy_xz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_0_xz[i] = -2.0 * g_y_xz_0_xz[i] * a_exp + 4.0 * g_xxy_xz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_0_yy[i] = -2.0 * g_y_xz_0_yy[i] * a_exp + 4.0 * g_xxy_xz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_0_yz[i] = -2.0 * g_y_xz_0_yz[i] * a_exp + 4.0 * g_xxy_xz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_0_zz[i] = -2.0 * g_y_xz_0_zz[i] * a_exp + 4.0 * g_xxy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xxy_yy_0_xx, g_xxy_yy_0_xy, g_xxy_yy_0_xz, g_xxy_yy_0_yy, g_xxy_yy_0_yz, g_xxy_yy_0_zz, g_xy_0_0_0_x_yy_0_xx, g_xy_0_0_0_x_yy_0_xy, g_xy_0_0_0_x_yy_0_xz, g_xy_0_0_0_x_yy_0_yy, g_xy_0_0_0_x_yy_0_yz, g_xy_0_0_0_x_yy_0_zz, g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_yy_0_xx[i] = -2.0 * g_y_yy_0_xx[i] * a_exp + 4.0 * g_xxy_yy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_0_xy[i] = -2.0 * g_y_yy_0_xy[i] * a_exp + 4.0 * g_xxy_yy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_0_xz[i] = -2.0 * g_y_yy_0_xz[i] * a_exp + 4.0 * g_xxy_yy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_0_yy[i] = -2.0 * g_y_yy_0_yy[i] * a_exp + 4.0 * g_xxy_yy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_0_yz[i] = -2.0 * g_y_yy_0_yz[i] * a_exp + 4.0 * g_xxy_yy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_0_zz[i] = -2.0 * g_y_yy_0_zz[i] * a_exp + 4.0 * g_xxy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_xxy_yz_0_xx, g_xxy_yz_0_xy, g_xxy_yz_0_xz, g_xxy_yz_0_yy, g_xxy_yz_0_yz, g_xxy_yz_0_zz, g_xy_0_0_0_x_yz_0_xx, g_xy_0_0_0_x_yz_0_xy, g_xy_0_0_0_x_yz_0_xz, g_xy_0_0_0_x_yz_0_yy, g_xy_0_0_0_x_yz_0_yz, g_xy_0_0_0_x_yz_0_zz, g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_yz_0_xx[i] = -2.0 * g_y_yz_0_xx[i] * a_exp + 4.0 * g_xxy_yz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_0_xy[i] = -2.0 * g_y_yz_0_xy[i] * a_exp + 4.0 * g_xxy_yz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_0_xz[i] = -2.0 * g_y_yz_0_xz[i] * a_exp + 4.0 * g_xxy_yz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_0_yy[i] = -2.0 * g_y_yz_0_yy[i] * a_exp + 4.0 * g_xxy_yz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_0_yz[i] = -2.0 * g_y_yz_0_yz[i] * a_exp + 4.0 * g_xxy_yz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_0_zz[i] = -2.0 * g_y_yz_0_zz[i] * a_exp + 4.0 * g_xxy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_xxy_zz_0_xx, g_xxy_zz_0_xy, g_xxy_zz_0_xz, g_xxy_zz_0_yy, g_xxy_zz_0_yz, g_xxy_zz_0_zz, g_xy_0_0_0_x_zz_0_xx, g_xy_0_0_0_x_zz_0_xy, g_xy_0_0_0_x_zz_0_xz, g_xy_0_0_0_x_zz_0_yy, g_xy_0_0_0_x_zz_0_yz, g_xy_0_0_0_x_zz_0_zz, g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_zz_0_xx[i] = -2.0 * g_y_zz_0_xx[i] * a_exp + 4.0 * g_xxy_zz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_0_xy[i] = -2.0 * g_y_zz_0_xy[i] * a_exp + 4.0 * g_xxy_zz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_0_xz[i] = -2.0 * g_y_zz_0_xz[i] * a_exp + 4.0 * g_xxy_zz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_0_yy[i] = -2.0 * g_y_zz_0_yy[i] * a_exp + 4.0 * g_xxy_zz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_0_yz[i] = -2.0 * g_y_zz_0_yz[i] * a_exp + 4.0 * g_xxy_zz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_0_zz[i] = -2.0 * g_y_zz_0_zz[i] * a_exp + 4.0 * g_xxy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_xy_0_0_0_y_xx_0_xx, g_xy_0_0_0_y_xx_0_xy, g_xy_0_0_0_y_xx_0_xz, g_xy_0_0_0_y_xx_0_yy, g_xy_0_0_0_y_xx_0_yz, g_xy_0_0_0_y_xx_0_zz, g_xyy_xx_0_xx, g_xyy_xx_0_xy, g_xyy_xx_0_xz, g_xyy_xx_0_yy, g_xyy_xx_0_yz, g_xyy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xx_0_xx[i] = -2.0 * g_x_xx_0_xx[i] * a_exp + 4.0 * g_xyy_xx_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_0_xy[i] = -2.0 * g_x_xx_0_xy[i] * a_exp + 4.0 * g_xyy_xx_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_0_xz[i] = -2.0 * g_x_xx_0_xz[i] * a_exp + 4.0 * g_xyy_xx_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_0_yy[i] = -2.0 * g_x_xx_0_yy[i] * a_exp + 4.0 * g_xyy_xx_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_0_yz[i] = -2.0 * g_x_xx_0_yz[i] * a_exp + 4.0 * g_xyy_xx_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_0_zz[i] = -2.0 * g_x_xx_0_zz[i] * a_exp + 4.0 * g_xyy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_xy_0_0_0_y_xy_0_xx, g_xy_0_0_0_y_xy_0_xy, g_xy_0_0_0_y_xy_0_xz, g_xy_0_0_0_y_xy_0_yy, g_xy_0_0_0_y_xy_0_yz, g_xy_0_0_0_y_xy_0_zz, g_xyy_xy_0_xx, g_xyy_xy_0_xy, g_xyy_xy_0_xz, g_xyy_xy_0_yy, g_xyy_xy_0_yz, g_xyy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xy_0_xx[i] = -2.0 * g_x_xy_0_xx[i] * a_exp + 4.0 * g_xyy_xy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_0_xy[i] = -2.0 * g_x_xy_0_xy[i] * a_exp + 4.0 * g_xyy_xy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_0_xz[i] = -2.0 * g_x_xy_0_xz[i] * a_exp + 4.0 * g_xyy_xy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_0_yy[i] = -2.0 * g_x_xy_0_yy[i] * a_exp + 4.0 * g_xyy_xy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_0_yz[i] = -2.0 * g_x_xy_0_yz[i] * a_exp + 4.0 * g_xyy_xy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_0_zz[i] = -2.0 * g_x_xy_0_zz[i] * a_exp + 4.0 * g_xyy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_xy_0_0_0_y_xz_0_xx, g_xy_0_0_0_y_xz_0_xy, g_xy_0_0_0_y_xz_0_xz, g_xy_0_0_0_y_xz_0_yy, g_xy_0_0_0_y_xz_0_yz, g_xy_0_0_0_y_xz_0_zz, g_xyy_xz_0_xx, g_xyy_xz_0_xy, g_xyy_xz_0_xz, g_xyy_xz_0_yy, g_xyy_xz_0_yz, g_xyy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xz_0_xx[i] = -2.0 * g_x_xz_0_xx[i] * a_exp + 4.0 * g_xyy_xz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_0_xy[i] = -2.0 * g_x_xz_0_xy[i] * a_exp + 4.0 * g_xyy_xz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_0_xz[i] = -2.0 * g_x_xz_0_xz[i] * a_exp + 4.0 * g_xyy_xz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_0_yy[i] = -2.0 * g_x_xz_0_yy[i] * a_exp + 4.0 * g_xyy_xz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_0_yz[i] = -2.0 * g_x_xz_0_yz[i] * a_exp + 4.0 * g_xyy_xz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_0_zz[i] = -2.0 * g_x_xz_0_zz[i] * a_exp + 4.0 * g_xyy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_xy_0_0_0_y_yy_0_xx, g_xy_0_0_0_y_yy_0_xy, g_xy_0_0_0_y_yy_0_xz, g_xy_0_0_0_y_yy_0_yy, g_xy_0_0_0_y_yy_0_yz, g_xy_0_0_0_y_yy_0_zz, g_xyy_yy_0_xx, g_xyy_yy_0_xy, g_xyy_yy_0_xz, g_xyy_yy_0_yy, g_xyy_yy_0_yz, g_xyy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_yy_0_xx[i] = -2.0 * g_x_yy_0_xx[i] * a_exp + 4.0 * g_xyy_yy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_0_xy[i] = -2.0 * g_x_yy_0_xy[i] * a_exp + 4.0 * g_xyy_yy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_0_xz[i] = -2.0 * g_x_yy_0_xz[i] * a_exp + 4.0 * g_xyy_yy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_0_yy[i] = -2.0 * g_x_yy_0_yy[i] * a_exp + 4.0 * g_xyy_yy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_0_yz[i] = -2.0 * g_x_yy_0_yz[i] * a_exp + 4.0 * g_xyy_yy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_0_zz[i] = -2.0 * g_x_yy_0_zz[i] * a_exp + 4.0 * g_xyy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_xy_0_0_0_y_yz_0_xx, g_xy_0_0_0_y_yz_0_xy, g_xy_0_0_0_y_yz_0_xz, g_xy_0_0_0_y_yz_0_yy, g_xy_0_0_0_y_yz_0_yz, g_xy_0_0_0_y_yz_0_zz, g_xyy_yz_0_xx, g_xyy_yz_0_xy, g_xyy_yz_0_xz, g_xyy_yz_0_yy, g_xyy_yz_0_yz, g_xyy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_yz_0_xx[i] = -2.0 * g_x_yz_0_xx[i] * a_exp + 4.0 * g_xyy_yz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_0_xy[i] = -2.0 * g_x_yz_0_xy[i] * a_exp + 4.0 * g_xyy_yz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_0_xz[i] = -2.0 * g_x_yz_0_xz[i] * a_exp + 4.0 * g_xyy_yz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_0_yy[i] = -2.0 * g_x_yz_0_yy[i] * a_exp + 4.0 * g_xyy_yz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_0_yz[i] = -2.0 * g_x_yz_0_yz[i] * a_exp + 4.0 * g_xyy_yz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_0_zz[i] = -2.0 * g_x_yz_0_zz[i] * a_exp + 4.0 * g_xyy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_xy_0_0_0_y_zz_0_xx, g_xy_0_0_0_y_zz_0_xy, g_xy_0_0_0_y_zz_0_xz, g_xy_0_0_0_y_zz_0_yy, g_xy_0_0_0_y_zz_0_yz, g_xy_0_0_0_y_zz_0_zz, g_xyy_zz_0_xx, g_xyy_zz_0_xy, g_xyy_zz_0_xz, g_xyy_zz_0_yy, g_xyy_zz_0_yz, g_xyy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_zz_0_xx[i] = -2.0 * g_x_zz_0_xx[i] * a_exp + 4.0 * g_xyy_zz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_0_xy[i] = -2.0 * g_x_zz_0_xy[i] * a_exp + 4.0 * g_xyy_zz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_0_xz[i] = -2.0 * g_x_zz_0_xz[i] * a_exp + 4.0 * g_xyy_zz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_0_yy[i] = -2.0 * g_x_zz_0_yy[i] * a_exp + 4.0 * g_xyy_zz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_0_yz[i] = -2.0 * g_x_zz_0_yz[i] * a_exp + 4.0 * g_xyy_zz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_0_zz[i] = -2.0 * g_x_zz_0_zz[i] * a_exp + 4.0 * g_xyy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xy_0_0_0_z_xx_0_xx, g_xy_0_0_0_z_xx_0_xy, g_xy_0_0_0_z_xx_0_xz, g_xy_0_0_0_z_xx_0_yy, g_xy_0_0_0_z_xx_0_yz, g_xy_0_0_0_z_xx_0_zz, g_xyz_xx_0_xx, g_xyz_xx_0_xy, g_xyz_xx_0_xz, g_xyz_xx_0_yy, g_xyz_xx_0_yz, g_xyz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xx_0_xx[i] = 4.0 * g_xyz_xx_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_0_xy[i] = 4.0 * g_xyz_xx_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_0_xz[i] = 4.0 * g_xyz_xx_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_0_yy[i] = 4.0 * g_xyz_xx_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_0_yz[i] = 4.0 * g_xyz_xx_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_0_zz[i] = 4.0 * g_xyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_xy_0_0_0_z_xy_0_xx, g_xy_0_0_0_z_xy_0_xy, g_xy_0_0_0_z_xy_0_xz, g_xy_0_0_0_z_xy_0_yy, g_xy_0_0_0_z_xy_0_yz, g_xy_0_0_0_z_xy_0_zz, g_xyz_xy_0_xx, g_xyz_xy_0_xy, g_xyz_xy_0_xz, g_xyz_xy_0_yy, g_xyz_xy_0_yz, g_xyz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xy_0_xx[i] = 4.0 * g_xyz_xy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_0_xy[i] = 4.0 * g_xyz_xy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_0_xz[i] = 4.0 * g_xyz_xy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_0_yy[i] = 4.0 * g_xyz_xy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_0_yz[i] = 4.0 * g_xyz_xy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_0_zz[i] = 4.0 * g_xyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xy_0_0_0_z_xz_0_xx, g_xy_0_0_0_z_xz_0_xy, g_xy_0_0_0_z_xz_0_xz, g_xy_0_0_0_z_xz_0_yy, g_xy_0_0_0_z_xz_0_yz, g_xy_0_0_0_z_xz_0_zz, g_xyz_xz_0_xx, g_xyz_xz_0_xy, g_xyz_xz_0_xz, g_xyz_xz_0_yy, g_xyz_xz_0_yz, g_xyz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xz_0_xx[i] = 4.0 * g_xyz_xz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_0_xy[i] = 4.0 * g_xyz_xz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_0_xz[i] = 4.0 * g_xyz_xz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_0_yy[i] = 4.0 * g_xyz_xz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_0_yz[i] = 4.0 * g_xyz_xz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_0_zz[i] = 4.0 * g_xyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_xy_0_0_0_z_yy_0_xx, g_xy_0_0_0_z_yy_0_xy, g_xy_0_0_0_z_yy_0_xz, g_xy_0_0_0_z_yy_0_yy, g_xy_0_0_0_z_yy_0_yz, g_xy_0_0_0_z_yy_0_zz, g_xyz_yy_0_xx, g_xyz_yy_0_xy, g_xyz_yy_0_xz, g_xyz_yy_0_yy, g_xyz_yy_0_yz, g_xyz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_yy_0_xx[i] = 4.0 * g_xyz_yy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_0_xy[i] = 4.0 * g_xyz_yy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_0_xz[i] = 4.0 * g_xyz_yy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_0_yy[i] = 4.0 * g_xyz_yy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_0_yz[i] = 4.0 * g_xyz_yy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_0_zz[i] = 4.0 * g_xyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_xy_0_0_0_z_yz_0_xx, g_xy_0_0_0_z_yz_0_xy, g_xy_0_0_0_z_yz_0_xz, g_xy_0_0_0_z_yz_0_yy, g_xy_0_0_0_z_yz_0_yz, g_xy_0_0_0_z_yz_0_zz, g_xyz_yz_0_xx, g_xyz_yz_0_xy, g_xyz_yz_0_xz, g_xyz_yz_0_yy, g_xyz_yz_0_yz, g_xyz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_yz_0_xx[i] = 4.0 * g_xyz_yz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_0_xy[i] = 4.0 * g_xyz_yz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_0_xz[i] = 4.0 * g_xyz_yz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_0_yy[i] = 4.0 * g_xyz_yz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_0_yz[i] = 4.0 * g_xyz_yz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_0_zz[i] = 4.0 * g_xyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_xy_0_0_0_z_zz_0_xx, g_xy_0_0_0_z_zz_0_xy, g_xy_0_0_0_z_zz_0_xz, g_xy_0_0_0_z_zz_0_yy, g_xy_0_0_0_z_zz_0_yz, g_xy_0_0_0_z_zz_0_zz, g_xyz_zz_0_xx, g_xyz_zz_0_xy, g_xyz_zz_0_xz, g_xyz_zz_0_yy, g_xyz_zz_0_yz, g_xyz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_zz_0_xx[i] = 4.0 * g_xyz_zz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_0_xy[i] = 4.0 * g_xyz_zz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_0_xz[i] = 4.0 * g_xyz_zz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_0_yy[i] = 4.0 * g_xyz_zz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_0_yz[i] = 4.0 * g_xyz_zz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_0_zz[i] = 4.0 * g_xyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xxz_xx_0_xx, g_xxz_xx_0_xy, g_xxz_xx_0_xz, g_xxz_xx_0_yy, g_xxz_xx_0_yz, g_xxz_xx_0_zz, g_xz_0_0_0_x_xx_0_xx, g_xz_0_0_0_x_xx_0_xy, g_xz_0_0_0_x_xx_0_xz, g_xz_0_0_0_x_xx_0_yy, g_xz_0_0_0_x_xx_0_yz, g_xz_0_0_0_x_xx_0_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xx_0_xx[i] = -2.0 * g_z_xx_0_xx[i] * a_exp + 4.0 * g_xxz_xx_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_0_xy[i] = -2.0 * g_z_xx_0_xy[i] * a_exp + 4.0 * g_xxz_xx_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_0_xz[i] = -2.0 * g_z_xx_0_xz[i] * a_exp + 4.0 * g_xxz_xx_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_0_yy[i] = -2.0 * g_z_xx_0_yy[i] * a_exp + 4.0 * g_xxz_xx_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_0_yz[i] = -2.0 * g_z_xx_0_yz[i] * a_exp + 4.0 * g_xxz_xx_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_0_zz[i] = -2.0 * g_z_xx_0_zz[i] * a_exp + 4.0 * g_xxz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xxz_xy_0_xx, g_xxz_xy_0_xy, g_xxz_xy_0_xz, g_xxz_xy_0_yy, g_xxz_xy_0_yz, g_xxz_xy_0_zz, g_xz_0_0_0_x_xy_0_xx, g_xz_0_0_0_x_xy_0_xy, g_xz_0_0_0_x_xy_0_xz, g_xz_0_0_0_x_xy_0_yy, g_xz_0_0_0_x_xy_0_yz, g_xz_0_0_0_x_xy_0_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xy_0_xx[i] = -2.0 * g_z_xy_0_xx[i] * a_exp + 4.0 * g_xxz_xy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_0_xy[i] = -2.0 * g_z_xy_0_xy[i] * a_exp + 4.0 * g_xxz_xy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_0_xz[i] = -2.0 * g_z_xy_0_xz[i] * a_exp + 4.0 * g_xxz_xy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_0_yy[i] = -2.0 * g_z_xy_0_yy[i] * a_exp + 4.0 * g_xxz_xy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_0_yz[i] = -2.0 * g_z_xy_0_yz[i] * a_exp + 4.0 * g_xxz_xy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_0_zz[i] = -2.0 * g_z_xy_0_zz[i] * a_exp + 4.0 * g_xxz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xxz_xz_0_xx, g_xxz_xz_0_xy, g_xxz_xz_0_xz, g_xxz_xz_0_yy, g_xxz_xz_0_yz, g_xxz_xz_0_zz, g_xz_0_0_0_x_xz_0_xx, g_xz_0_0_0_x_xz_0_xy, g_xz_0_0_0_x_xz_0_xz, g_xz_0_0_0_x_xz_0_yy, g_xz_0_0_0_x_xz_0_yz, g_xz_0_0_0_x_xz_0_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xz_0_xx[i] = -2.0 * g_z_xz_0_xx[i] * a_exp + 4.0 * g_xxz_xz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_0_xy[i] = -2.0 * g_z_xz_0_xy[i] * a_exp + 4.0 * g_xxz_xz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_0_xz[i] = -2.0 * g_z_xz_0_xz[i] * a_exp + 4.0 * g_xxz_xz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_0_yy[i] = -2.0 * g_z_xz_0_yy[i] * a_exp + 4.0 * g_xxz_xz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_0_yz[i] = -2.0 * g_z_xz_0_yz[i] * a_exp + 4.0 * g_xxz_xz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_0_zz[i] = -2.0 * g_z_xz_0_zz[i] * a_exp + 4.0 * g_xxz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_xxz_yy_0_xx, g_xxz_yy_0_xy, g_xxz_yy_0_xz, g_xxz_yy_0_yy, g_xxz_yy_0_yz, g_xxz_yy_0_zz, g_xz_0_0_0_x_yy_0_xx, g_xz_0_0_0_x_yy_0_xy, g_xz_0_0_0_x_yy_0_xz, g_xz_0_0_0_x_yy_0_yy, g_xz_0_0_0_x_yy_0_yz, g_xz_0_0_0_x_yy_0_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_yy_0_xx[i] = -2.0 * g_z_yy_0_xx[i] * a_exp + 4.0 * g_xxz_yy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_0_xy[i] = -2.0 * g_z_yy_0_xy[i] * a_exp + 4.0 * g_xxz_yy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_0_xz[i] = -2.0 * g_z_yy_0_xz[i] * a_exp + 4.0 * g_xxz_yy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_0_yy[i] = -2.0 * g_z_yy_0_yy[i] * a_exp + 4.0 * g_xxz_yy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_0_yz[i] = -2.0 * g_z_yy_0_yz[i] * a_exp + 4.0 * g_xxz_yy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_0_zz[i] = -2.0 * g_z_yy_0_zz[i] * a_exp + 4.0 * g_xxz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_xxz_yz_0_xx, g_xxz_yz_0_xy, g_xxz_yz_0_xz, g_xxz_yz_0_yy, g_xxz_yz_0_yz, g_xxz_yz_0_zz, g_xz_0_0_0_x_yz_0_xx, g_xz_0_0_0_x_yz_0_xy, g_xz_0_0_0_x_yz_0_xz, g_xz_0_0_0_x_yz_0_yy, g_xz_0_0_0_x_yz_0_yz, g_xz_0_0_0_x_yz_0_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_yz_0_xx[i] = -2.0 * g_z_yz_0_xx[i] * a_exp + 4.0 * g_xxz_yz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_0_xy[i] = -2.0 * g_z_yz_0_xy[i] * a_exp + 4.0 * g_xxz_yz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_0_xz[i] = -2.0 * g_z_yz_0_xz[i] * a_exp + 4.0 * g_xxz_yz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_0_yy[i] = -2.0 * g_z_yz_0_yy[i] * a_exp + 4.0 * g_xxz_yz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_0_yz[i] = -2.0 * g_z_yz_0_yz[i] * a_exp + 4.0 * g_xxz_yz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_0_zz[i] = -2.0 * g_z_yz_0_zz[i] * a_exp + 4.0 * g_xxz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_xxz_zz_0_xx, g_xxz_zz_0_xy, g_xxz_zz_0_xz, g_xxz_zz_0_yy, g_xxz_zz_0_yz, g_xxz_zz_0_zz, g_xz_0_0_0_x_zz_0_xx, g_xz_0_0_0_x_zz_0_xy, g_xz_0_0_0_x_zz_0_xz, g_xz_0_0_0_x_zz_0_yy, g_xz_0_0_0_x_zz_0_yz, g_xz_0_0_0_x_zz_0_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_zz_0_xx[i] = -2.0 * g_z_zz_0_xx[i] * a_exp + 4.0 * g_xxz_zz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_0_xy[i] = -2.0 * g_z_zz_0_xy[i] * a_exp + 4.0 * g_xxz_zz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_0_xz[i] = -2.0 * g_z_zz_0_xz[i] * a_exp + 4.0 * g_xxz_zz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_0_yy[i] = -2.0 * g_z_zz_0_yy[i] * a_exp + 4.0 * g_xxz_zz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_0_yz[i] = -2.0 * g_z_zz_0_yz[i] * a_exp + 4.0 * g_xxz_zz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_0_zz[i] = -2.0 * g_z_zz_0_zz[i] * a_exp + 4.0 * g_xxz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_xyz_xx_0_xx, g_xyz_xx_0_xy, g_xyz_xx_0_xz, g_xyz_xx_0_yy, g_xyz_xx_0_yz, g_xyz_xx_0_zz, g_xz_0_0_0_y_xx_0_xx, g_xz_0_0_0_y_xx_0_xy, g_xz_0_0_0_y_xx_0_xz, g_xz_0_0_0_y_xx_0_yy, g_xz_0_0_0_y_xx_0_yz, g_xz_0_0_0_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xx_0_xx[i] = 4.0 * g_xyz_xx_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_0_xy[i] = 4.0 * g_xyz_xx_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_0_xz[i] = 4.0 * g_xyz_xx_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_0_yy[i] = 4.0 * g_xyz_xx_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_0_yz[i] = 4.0 * g_xyz_xx_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_0_zz[i] = 4.0 * g_xyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_xyz_xy_0_xx, g_xyz_xy_0_xy, g_xyz_xy_0_xz, g_xyz_xy_0_yy, g_xyz_xy_0_yz, g_xyz_xy_0_zz, g_xz_0_0_0_y_xy_0_xx, g_xz_0_0_0_y_xy_0_xy, g_xz_0_0_0_y_xy_0_xz, g_xz_0_0_0_y_xy_0_yy, g_xz_0_0_0_y_xy_0_yz, g_xz_0_0_0_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xy_0_xx[i] = 4.0 * g_xyz_xy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_0_xy[i] = 4.0 * g_xyz_xy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_0_xz[i] = 4.0 * g_xyz_xy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_0_yy[i] = 4.0 * g_xyz_xy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_0_yz[i] = 4.0 * g_xyz_xy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_0_zz[i] = 4.0 * g_xyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_xyz_xz_0_xx, g_xyz_xz_0_xy, g_xyz_xz_0_xz, g_xyz_xz_0_yy, g_xyz_xz_0_yz, g_xyz_xz_0_zz, g_xz_0_0_0_y_xz_0_xx, g_xz_0_0_0_y_xz_0_xy, g_xz_0_0_0_y_xz_0_xz, g_xz_0_0_0_y_xz_0_yy, g_xz_0_0_0_y_xz_0_yz, g_xz_0_0_0_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xz_0_xx[i] = 4.0 * g_xyz_xz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_0_xy[i] = 4.0 * g_xyz_xz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_0_xz[i] = 4.0 * g_xyz_xz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_0_yy[i] = 4.0 * g_xyz_xz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_0_yz[i] = 4.0 * g_xyz_xz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_0_zz[i] = 4.0 * g_xyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_xyz_yy_0_xx, g_xyz_yy_0_xy, g_xyz_yy_0_xz, g_xyz_yy_0_yy, g_xyz_yy_0_yz, g_xyz_yy_0_zz, g_xz_0_0_0_y_yy_0_xx, g_xz_0_0_0_y_yy_0_xy, g_xz_0_0_0_y_yy_0_xz, g_xz_0_0_0_y_yy_0_yy, g_xz_0_0_0_y_yy_0_yz, g_xz_0_0_0_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_yy_0_xx[i] = 4.0 * g_xyz_yy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_0_xy[i] = 4.0 * g_xyz_yy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_0_xz[i] = 4.0 * g_xyz_yy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_0_yy[i] = 4.0 * g_xyz_yy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_0_yz[i] = 4.0 * g_xyz_yy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_0_zz[i] = 4.0 * g_xyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_xyz_yz_0_xx, g_xyz_yz_0_xy, g_xyz_yz_0_xz, g_xyz_yz_0_yy, g_xyz_yz_0_yz, g_xyz_yz_0_zz, g_xz_0_0_0_y_yz_0_xx, g_xz_0_0_0_y_yz_0_xy, g_xz_0_0_0_y_yz_0_xz, g_xz_0_0_0_y_yz_0_yy, g_xz_0_0_0_y_yz_0_yz, g_xz_0_0_0_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_yz_0_xx[i] = 4.0 * g_xyz_yz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_0_xy[i] = 4.0 * g_xyz_yz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_0_xz[i] = 4.0 * g_xyz_yz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_0_yy[i] = 4.0 * g_xyz_yz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_0_yz[i] = 4.0 * g_xyz_yz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_0_zz[i] = 4.0 * g_xyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_xyz_zz_0_xx, g_xyz_zz_0_xy, g_xyz_zz_0_xz, g_xyz_zz_0_yy, g_xyz_zz_0_yz, g_xyz_zz_0_zz, g_xz_0_0_0_y_zz_0_xx, g_xz_0_0_0_y_zz_0_xy, g_xz_0_0_0_y_zz_0_xz, g_xz_0_0_0_y_zz_0_yy, g_xz_0_0_0_y_zz_0_yz, g_xz_0_0_0_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_zz_0_xx[i] = 4.0 * g_xyz_zz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_0_xy[i] = 4.0 * g_xyz_zz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_0_xz[i] = 4.0 * g_xyz_zz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_0_yy[i] = 4.0 * g_xyz_zz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_0_yz[i] = 4.0 * g_xyz_zz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_0_zz[i] = 4.0 * g_xyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_xz_0_0_0_z_xx_0_xx, g_xz_0_0_0_z_xx_0_xy, g_xz_0_0_0_z_xx_0_xz, g_xz_0_0_0_z_xx_0_yy, g_xz_0_0_0_z_xx_0_yz, g_xz_0_0_0_z_xx_0_zz, g_xzz_xx_0_xx, g_xzz_xx_0_xy, g_xzz_xx_0_xz, g_xzz_xx_0_yy, g_xzz_xx_0_yz, g_xzz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xx_0_xx[i] = -2.0 * g_x_xx_0_xx[i] * a_exp + 4.0 * g_xzz_xx_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_0_xy[i] = -2.0 * g_x_xx_0_xy[i] * a_exp + 4.0 * g_xzz_xx_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_0_xz[i] = -2.0 * g_x_xx_0_xz[i] * a_exp + 4.0 * g_xzz_xx_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_0_yy[i] = -2.0 * g_x_xx_0_yy[i] * a_exp + 4.0 * g_xzz_xx_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_0_yz[i] = -2.0 * g_x_xx_0_yz[i] * a_exp + 4.0 * g_xzz_xx_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_0_zz[i] = -2.0 * g_x_xx_0_zz[i] * a_exp + 4.0 * g_xzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_xz_0_0_0_z_xy_0_xx, g_xz_0_0_0_z_xy_0_xy, g_xz_0_0_0_z_xy_0_xz, g_xz_0_0_0_z_xy_0_yy, g_xz_0_0_0_z_xy_0_yz, g_xz_0_0_0_z_xy_0_zz, g_xzz_xy_0_xx, g_xzz_xy_0_xy, g_xzz_xy_0_xz, g_xzz_xy_0_yy, g_xzz_xy_0_yz, g_xzz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xy_0_xx[i] = -2.0 * g_x_xy_0_xx[i] * a_exp + 4.0 * g_xzz_xy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_0_xy[i] = -2.0 * g_x_xy_0_xy[i] * a_exp + 4.0 * g_xzz_xy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_0_xz[i] = -2.0 * g_x_xy_0_xz[i] * a_exp + 4.0 * g_xzz_xy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_0_yy[i] = -2.0 * g_x_xy_0_yy[i] * a_exp + 4.0 * g_xzz_xy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_0_yz[i] = -2.0 * g_x_xy_0_yz[i] * a_exp + 4.0 * g_xzz_xy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_0_zz[i] = -2.0 * g_x_xy_0_zz[i] * a_exp + 4.0 * g_xzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_xz_0_0_0_z_xz_0_xx, g_xz_0_0_0_z_xz_0_xy, g_xz_0_0_0_z_xz_0_xz, g_xz_0_0_0_z_xz_0_yy, g_xz_0_0_0_z_xz_0_yz, g_xz_0_0_0_z_xz_0_zz, g_xzz_xz_0_xx, g_xzz_xz_0_xy, g_xzz_xz_0_xz, g_xzz_xz_0_yy, g_xzz_xz_0_yz, g_xzz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xz_0_xx[i] = -2.0 * g_x_xz_0_xx[i] * a_exp + 4.0 * g_xzz_xz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_0_xy[i] = -2.0 * g_x_xz_0_xy[i] * a_exp + 4.0 * g_xzz_xz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_0_xz[i] = -2.0 * g_x_xz_0_xz[i] * a_exp + 4.0 * g_xzz_xz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_0_yy[i] = -2.0 * g_x_xz_0_yy[i] * a_exp + 4.0 * g_xzz_xz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_0_yz[i] = -2.0 * g_x_xz_0_yz[i] * a_exp + 4.0 * g_xzz_xz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_0_zz[i] = -2.0 * g_x_xz_0_zz[i] * a_exp + 4.0 * g_xzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_xz_0_0_0_z_yy_0_xx, g_xz_0_0_0_z_yy_0_xy, g_xz_0_0_0_z_yy_0_xz, g_xz_0_0_0_z_yy_0_yy, g_xz_0_0_0_z_yy_0_yz, g_xz_0_0_0_z_yy_0_zz, g_xzz_yy_0_xx, g_xzz_yy_0_xy, g_xzz_yy_0_xz, g_xzz_yy_0_yy, g_xzz_yy_0_yz, g_xzz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_yy_0_xx[i] = -2.0 * g_x_yy_0_xx[i] * a_exp + 4.0 * g_xzz_yy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_0_xy[i] = -2.0 * g_x_yy_0_xy[i] * a_exp + 4.0 * g_xzz_yy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_0_xz[i] = -2.0 * g_x_yy_0_xz[i] * a_exp + 4.0 * g_xzz_yy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_0_yy[i] = -2.0 * g_x_yy_0_yy[i] * a_exp + 4.0 * g_xzz_yy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_0_yz[i] = -2.0 * g_x_yy_0_yz[i] * a_exp + 4.0 * g_xzz_yy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_0_zz[i] = -2.0 * g_x_yy_0_zz[i] * a_exp + 4.0 * g_xzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_xz_0_0_0_z_yz_0_xx, g_xz_0_0_0_z_yz_0_xy, g_xz_0_0_0_z_yz_0_xz, g_xz_0_0_0_z_yz_0_yy, g_xz_0_0_0_z_yz_0_yz, g_xz_0_0_0_z_yz_0_zz, g_xzz_yz_0_xx, g_xzz_yz_0_xy, g_xzz_yz_0_xz, g_xzz_yz_0_yy, g_xzz_yz_0_yz, g_xzz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_yz_0_xx[i] = -2.0 * g_x_yz_0_xx[i] * a_exp + 4.0 * g_xzz_yz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_0_xy[i] = -2.0 * g_x_yz_0_xy[i] * a_exp + 4.0 * g_xzz_yz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_0_xz[i] = -2.0 * g_x_yz_0_xz[i] * a_exp + 4.0 * g_xzz_yz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_0_yy[i] = -2.0 * g_x_yz_0_yy[i] * a_exp + 4.0 * g_xzz_yz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_0_yz[i] = -2.0 * g_x_yz_0_yz[i] * a_exp + 4.0 * g_xzz_yz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_0_zz[i] = -2.0 * g_x_yz_0_zz[i] * a_exp + 4.0 * g_xzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_xz_0_0_0_z_zz_0_xx, g_xz_0_0_0_z_zz_0_xy, g_xz_0_0_0_z_zz_0_xz, g_xz_0_0_0_z_zz_0_yy, g_xz_0_0_0_z_zz_0_yz, g_xz_0_0_0_z_zz_0_zz, g_xzz_zz_0_xx, g_xzz_zz_0_xy, g_xzz_zz_0_xz, g_xzz_zz_0_yy, g_xzz_zz_0_yz, g_xzz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_zz_0_xx[i] = -2.0 * g_x_zz_0_xx[i] * a_exp + 4.0 * g_xzz_zz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_0_xy[i] = -2.0 * g_x_zz_0_xy[i] * a_exp + 4.0 * g_xzz_zz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_0_xz[i] = -2.0 * g_x_zz_0_xz[i] * a_exp + 4.0 * g_xzz_zz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_0_yy[i] = -2.0 * g_x_zz_0_yy[i] * a_exp + 4.0 * g_xzz_zz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_0_yz[i] = -2.0 * g_x_zz_0_yz[i] * a_exp + 4.0 * g_xzz_zz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_0_zz[i] = -2.0 * g_x_zz_0_zz[i] * a_exp + 4.0 * g_xzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_xyy_xx_0_xx, g_xyy_xx_0_xy, g_xyy_xx_0_xz, g_xyy_xx_0_yy, g_xyy_xx_0_yz, g_xyy_xx_0_zz, g_yy_0_0_0_x_xx_0_xx, g_yy_0_0_0_x_xx_0_xy, g_yy_0_0_0_x_xx_0_xz, g_yy_0_0_0_x_xx_0_yy, g_yy_0_0_0_x_xx_0_yz, g_yy_0_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xx_0_xx[i] = -2.0 * g_x_xx_0_xx[i] * a_exp + 4.0 * g_xyy_xx_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_0_xy[i] = -2.0 * g_x_xx_0_xy[i] * a_exp + 4.0 * g_xyy_xx_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_0_xz[i] = -2.0 * g_x_xx_0_xz[i] * a_exp + 4.0 * g_xyy_xx_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_0_yy[i] = -2.0 * g_x_xx_0_yy[i] * a_exp + 4.0 * g_xyy_xx_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_0_yz[i] = -2.0 * g_x_xx_0_yz[i] * a_exp + 4.0 * g_xyy_xx_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_0_zz[i] = -2.0 * g_x_xx_0_zz[i] * a_exp + 4.0 * g_xyy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_xyy_xy_0_xx, g_xyy_xy_0_xy, g_xyy_xy_0_xz, g_xyy_xy_0_yy, g_xyy_xy_0_yz, g_xyy_xy_0_zz, g_yy_0_0_0_x_xy_0_xx, g_yy_0_0_0_x_xy_0_xy, g_yy_0_0_0_x_xy_0_xz, g_yy_0_0_0_x_xy_0_yy, g_yy_0_0_0_x_xy_0_yz, g_yy_0_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xy_0_xx[i] = -2.0 * g_x_xy_0_xx[i] * a_exp + 4.0 * g_xyy_xy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_0_xy[i] = -2.0 * g_x_xy_0_xy[i] * a_exp + 4.0 * g_xyy_xy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_0_xz[i] = -2.0 * g_x_xy_0_xz[i] * a_exp + 4.0 * g_xyy_xy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_0_yy[i] = -2.0 * g_x_xy_0_yy[i] * a_exp + 4.0 * g_xyy_xy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_0_yz[i] = -2.0 * g_x_xy_0_yz[i] * a_exp + 4.0 * g_xyy_xy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_0_zz[i] = -2.0 * g_x_xy_0_zz[i] * a_exp + 4.0 * g_xyy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_xyy_xz_0_xx, g_xyy_xz_0_xy, g_xyy_xz_0_xz, g_xyy_xz_0_yy, g_xyy_xz_0_yz, g_xyy_xz_0_zz, g_yy_0_0_0_x_xz_0_xx, g_yy_0_0_0_x_xz_0_xy, g_yy_0_0_0_x_xz_0_xz, g_yy_0_0_0_x_xz_0_yy, g_yy_0_0_0_x_xz_0_yz, g_yy_0_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xz_0_xx[i] = -2.0 * g_x_xz_0_xx[i] * a_exp + 4.0 * g_xyy_xz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_0_xy[i] = -2.0 * g_x_xz_0_xy[i] * a_exp + 4.0 * g_xyy_xz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_0_xz[i] = -2.0 * g_x_xz_0_xz[i] * a_exp + 4.0 * g_xyy_xz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_0_yy[i] = -2.0 * g_x_xz_0_yy[i] * a_exp + 4.0 * g_xyy_xz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_0_yz[i] = -2.0 * g_x_xz_0_yz[i] * a_exp + 4.0 * g_xyy_xz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_0_zz[i] = -2.0 * g_x_xz_0_zz[i] * a_exp + 4.0 * g_xyy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_xyy_yy_0_xx, g_xyy_yy_0_xy, g_xyy_yy_0_xz, g_xyy_yy_0_yy, g_xyy_yy_0_yz, g_xyy_yy_0_zz, g_yy_0_0_0_x_yy_0_xx, g_yy_0_0_0_x_yy_0_xy, g_yy_0_0_0_x_yy_0_xz, g_yy_0_0_0_x_yy_0_yy, g_yy_0_0_0_x_yy_0_yz, g_yy_0_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_yy_0_xx[i] = -2.0 * g_x_yy_0_xx[i] * a_exp + 4.0 * g_xyy_yy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_0_xy[i] = -2.0 * g_x_yy_0_xy[i] * a_exp + 4.0 * g_xyy_yy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_0_xz[i] = -2.0 * g_x_yy_0_xz[i] * a_exp + 4.0 * g_xyy_yy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_0_yy[i] = -2.0 * g_x_yy_0_yy[i] * a_exp + 4.0 * g_xyy_yy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_0_yz[i] = -2.0 * g_x_yy_0_yz[i] * a_exp + 4.0 * g_xyy_yy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_0_zz[i] = -2.0 * g_x_yy_0_zz[i] * a_exp + 4.0 * g_xyy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_xyy_yz_0_xx, g_xyy_yz_0_xy, g_xyy_yz_0_xz, g_xyy_yz_0_yy, g_xyy_yz_0_yz, g_xyy_yz_0_zz, g_yy_0_0_0_x_yz_0_xx, g_yy_0_0_0_x_yz_0_xy, g_yy_0_0_0_x_yz_0_xz, g_yy_0_0_0_x_yz_0_yy, g_yy_0_0_0_x_yz_0_yz, g_yy_0_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_yz_0_xx[i] = -2.0 * g_x_yz_0_xx[i] * a_exp + 4.0 * g_xyy_yz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_0_xy[i] = -2.0 * g_x_yz_0_xy[i] * a_exp + 4.0 * g_xyy_yz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_0_xz[i] = -2.0 * g_x_yz_0_xz[i] * a_exp + 4.0 * g_xyy_yz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_0_yy[i] = -2.0 * g_x_yz_0_yy[i] * a_exp + 4.0 * g_xyy_yz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_0_yz[i] = -2.0 * g_x_yz_0_yz[i] * a_exp + 4.0 * g_xyy_yz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_0_zz[i] = -2.0 * g_x_yz_0_zz[i] * a_exp + 4.0 * g_xyy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_xyy_zz_0_xx, g_xyy_zz_0_xy, g_xyy_zz_0_xz, g_xyy_zz_0_yy, g_xyy_zz_0_yz, g_xyy_zz_0_zz, g_yy_0_0_0_x_zz_0_xx, g_yy_0_0_0_x_zz_0_xy, g_yy_0_0_0_x_zz_0_xz, g_yy_0_0_0_x_zz_0_yy, g_yy_0_0_0_x_zz_0_yz, g_yy_0_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_zz_0_xx[i] = -2.0 * g_x_zz_0_xx[i] * a_exp + 4.0 * g_xyy_zz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_0_xy[i] = -2.0 * g_x_zz_0_xy[i] * a_exp + 4.0 * g_xyy_zz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_0_xz[i] = -2.0 * g_x_zz_0_xz[i] * a_exp + 4.0 * g_xyy_zz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_0_yy[i] = -2.0 * g_x_zz_0_yy[i] * a_exp + 4.0 * g_xyy_zz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_0_yz[i] = -2.0 * g_x_zz_0_yz[i] * a_exp + 4.0 * g_xyy_zz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_0_zz[i] = -2.0 * g_x_zz_0_zz[i] * a_exp + 4.0 * g_xyy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz, g_yy_0_0_0_y_xx_0_xx, g_yy_0_0_0_y_xx_0_xy, g_yy_0_0_0_y_xx_0_xz, g_yy_0_0_0_y_xx_0_yy, g_yy_0_0_0_y_xx_0_yz, g_yy_0_0_0_y_xx_0_zz, g_yyy_xx_0_xx, g_yyy_xx_0_xy, g_yyy_xx_0_xz, g_yyy_xx_0_yy, g_yyy_xx_0_yz, g_yyy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xx_0_xx[i] = -6.0 * g_y_xx_0_xx[i] * a_exp + 4.0 * g_yyy_xx_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_0_xy[i] = -6.0 * g_y_xx_0_xy[i] * a_exp + 4.0 * g_yyy_xx_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_0_xz[i] = -6.0 * g_y_xx_0_xz[i] * a_exp + 4.0 * g_yyy_xx_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_0_yy[i] = -6.0 * g_y_xx_0_yy[i] * a_exp + 4.0 * g_yyy_xx_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_0_yz[i] = -6.0 * g_y_xx_0_yz[i] * a_exp + 4.0 * g_yyy_xx_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_0_zz[i] = -6.0 * g_y_xx_0_zz[i] * a_exp + 4.0 * g_yyy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz, g_yy_0_0_0_y_xy_0_xx, g_yy_0_0_0_y_xy_0_xy, g_yy_0_0_0_y_xy_0_xz, g_yy_0_0_0_y_xy_0_yy, g_yy_0_0_0_y_xy_0_yz, g_yy_0_0_0_y_xy_0_zz, g_yyy_xy_0_xx, g_yyy_xy_0_xy, g_yyy_xy_0_xz, g_yyy_xy_0_yy, g_yyy_xy_0_yz, g_yyy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xy_0_xx[i] = -6.0 * g_y_xy_0_xx[i] * a_exp + 4.0 * g_yyy_xy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_0_xy[i] = -6.0 * g_y_xy_0_xy[i] * a_exp + 4.0 * g_yyy_xy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_0_xz[i] = -6.0 * g_y_xy_0_xz[i] * a_exp + 4.0 * g_yyy_xy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_0_yy[i] = -6.0 * g_y_xy_0_yy[i] * a_exp + 4.0 * g_yyy_xy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_0_yz[i] = -6.0 * g_y_xy_0_yz[i] * a_exp + 4.0 * g_yyy_xy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_0_zz[i] = -6.0 * g_y_xy_0_zz[i] * a_exp + 4.0 * g_yyy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz, g_yy_0_0_0_y_xz_0_xx, g_yy_0_0_0_y_xz_0_xy, g_yy_0_0_0_y_xz_0_xz, g_yy_0_0_0_y_xz_0_yy, g_yy_0_0_0_y_xz_0_yz, g_yy_0_0_0_y_xz_0_zz, g_yyy_xz_0_xx, g_yyy_xz_0_xy, g_yyy_xz_0_xz, g_yyy_xz_0_yy, g_yyy_xz_0_yz, g_yyy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xz_0_xx[i] = -6.0 * g_y_xz_0_xx[i] * a_exp + 4.0 * g_yyy_xz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_0_xy[i] = -6.0 * g_y_xz_0_xy[i] * a_exp + 4.0 * g_yyy_xz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_0_xz[i] = -6.0 * g_y_xz_0_xz[i] * a_exp + 4.0 * g_yyy_xz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_0_yy[i] = -6.0 * g_y_xz_0_yy[i] * a_exp + 4.0 * g_yyy_xz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_0_yz[i] = -6.0 * g_y_xz_0_yz[i] * a_exp + 4.0 * g_yyy_xz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_0_zz[i] = -6.0 * g_y_xz_0_zz[i] * a_exp + 4.0 * g_yyy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz, g_yy_0_0_0_y_yy_0_xx, g_yy_0_0_0_y_yy_0_xy, g_yy_0_0_0_y_yy_0_xz, g_yy_0_0_0_y_yy_0_yy, g_yy_0_0_0_y_yy_0_yz, g_yy_0_0_0_y_yy_0_zz, g_yyy_yy_0_xx, g_yyy_yy_0_xy, g_yyy_yy_0_xz, g_yyy_yy_0_yy, g_yyy_yy_0_yz, g_yyy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_yy_0_xx[i] = -6.0 * g_y_yy_0_xx[i] * a_exp + 4.0 * g_yyy_yy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_0_xy[i] = -6.0 * g_y_yy_0_xy[i] * a_exp + 4.0 * g_yyy_yy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_0_xz[i] = -6.0 * g_y_yy_0_xz[i] * a_exp + 4.0 * g_yyy_yy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_0_yy[i] = -6.0 * g_y_yy_0_yy[i] * a_exp + 4.0 * g_yyy_yy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_0_yz[i] = -6.0 * g_y_yy_0_yz[i] * a_exp + 4.0 * g_yyy_yy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_0_zz[i] = -6.0 * g_y_yy_0_zz[i] * a_exp + 4.0 * g_yyy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz, g_yy_0_0_0_y_yz_0_xx, g_yy_0_0_0_y_yz_0_xy, g_yy_0_0_0_y_yz_0_xz, g_yy_0_0_0_y_yz_0_yy, g_yy_0_0_0_y_yz_0_yz, g_yy_0_0_0_y_yz_0_zz, g_yyy_yz_0_xx, g_yyy_yz_0_xy, g_yyy_yz_0_xz, g_yyy_yz_0_yy, g_yyy_yz_0_yz, g_yyy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_yz_0_xx[i] = -6.0 * g_y_yz_0_xx[i] * a_exp + 4.0 * g_yyy_yz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_0_xy[i] = -6.0 * g_y_yz_0_xy[i] * a_exp + 4.0 * g_yyy_yz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_0_xz[i] = -6.0 * g_y_yz_0_xz[i] * a_exp + 4.0 * g_yyy_yz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_0_yy[i] = -6.0 * g_y_yz_0_yy[i] * a_exp + 4.0 * g_yyy_yz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_0_yz[i] = -6.0 * g_y_yz_0_yz[i] * a_exp + 4.0 * g_yyy_yz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_0_zz[i] = -6.0 * g_y_yz_0_zz[i] * a_exp + 4.0 * g_yyy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz, g_yy_0_0_0_y_zz_0_xx, g_yy_0_0_0_y_zz_0_xy, g_yy_0_0_0_y_zz_0_xz, g_yy_0_0_0_y_zz_0_yy, g_yy_0_0_0_y_zz_0_yz, g_yy_0_0_0_y_zz_0_zz, g_yyy_zz_0_xx, g_yyy_zz_0_xy, g_yyy_zz_0_xz, g_yyy_zz_0_yy, g_yyy_zz_0_yz, g_yyy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_zz_0_xx[i] = -6.0 * g_y_zz_0_xx[i] * a_exp + 4.0 * g_yyy_zz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_0_xy[i] = -6.0 * g_y_zz_0_xy[i] * a_exp + 4.0 * g_yyy_zz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_0_xz[i] = -6.0 * g_y_zz_0_xz[i] * a_exp + 4.0 * g_yyy_zz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_0_yy[i] = -6.0 * g_y_zz_0_yy[i] * a_exp + 4.0 * g_yyy_zz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_0_yz[i] = -6.0 * g_y_zz_0_yz[i] * a_exp + 4.0 * g_yyy_zz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_0_zz[i] = -6.0 * g_y_zz_0_zz[i] * a_exp + 4.0 * g_yyy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_yy_0_0_0_z_xx_0_xx, g_yy_0_0_0_z_xx_0_xy, g_yy_0_0_0_z_xx_0_xz, g_yy_0_0_0_z_xx_0_yy, g_yy_0_0_0_z_xx_0_yz, g_yy_0_0_0_z_xx_0_zz, g_yyz_xx_0_xx, g_yyz_xx_0_xy, g_yyz_xx_0_xz, g_yyz_xx_0_yy, g_yyz_xx_0_yz, g_yyz_xx_0_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xx_0_xx[i] = -2.0 * g_z_xx_0_xx[i] * a_exp + 4.0 * g_yyz_xx_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_0_xy[i] = -2.0 * g_z_xx_0_xy[i] * a_exp + 4.0 * g_yyz_xx_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_0_xz[i] = -2.0 * g_z_xx_0_xz[i] * a_exp + 4.0 * g_yyz_xx_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_0_yy[i] = -2.0 * g_z_xx_0_yy[i] * a_exp + 4.0 * g_yyz_xx_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_0_yz[i] = -2.0 * g_z_xx_0_yz[i] * a_exp + 4.0 * g_yyz_xx_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_0_zz[i] = -2.0 * g_z_xx_0_zz[i] * a_exp + 4.0 * g_yyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_yy_0_0_0_z_xy_0_xx, g_yy_0_0_0_z_xy_0_xy, g_yy_0_0_0_z_xy_0_xz, g_yy_0_0_0_z_xy_0_yy, g_yy_0_0_0_z_xy_0_yz, g_yy_0_0_0_z_xy_0_zz, g_yyz_xy_0_xx, g_yyz_xy_0_xy, g_yyz_xy_0_xz, g_yyz_xy_0_yy, g_yyz_xy_0_yz, g_yyz_xy_0_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xy_0_xx[i] = -2.0 * g_z_xy_0_xx[i] * a_exp + 4.0 * g_yyz_xy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_0_xy[i] = -2.0 * g_z_xy_0_xy[i] * a_exp + 4.0 * g_yyz_xy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_0_xz[i] = -2.0 * g_z_xy_0_xz[i] * a_exp + 4.0 * g_yyz_xy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_0_yy[i] = -2.0 * g_z_xy_0_yy[i] * a_exp + 4.0 * g_yyz_xy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_0_yz[i] = -2.0 * g_z_xy_0_yz[i] * a_exp + 4.0 * g_yyz_xy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_0_zz[i] = -2.0 * g_z_xy_0_zz[i] * a_exp + 4.0 * g_yyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_yy_0_0_0_z_xz_0_xx, g_yy_0_0_0_z_xz_0_xy, g_yy_0_0_0_z_xz_0_xz, g_yy_0_0_0_z_xz_0_yy, g_yy_0_0_0_z_xz_0_yz, g_yy_0_0_0_z_xz_0_zz, g_yyz_xz_0_xx, g_yyz_xz_0_xy, g_yyz_xz_0_xz, g_yyz_xz_0_yy, g_yyz_xz_0_yz, g_yyz_xz_0_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xz_0_xx[i] = -2.0 * g_z_xz_0_xx[i] * a_exp + 4.0 * g_yyz_xz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_0_xy[i] = -2.0 * g_z_xz_0_xy[i] * a_exp + 4.0 * g_yyz_xz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_0_xz[i] = -2.0 * g_z_xz_0_xz[i] * a_exp + 4.0 * g_yyz_xz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_0_yy[i] = -2.0 * g_z_xz_0_yy[i] * a_exp + 4.0 * g_yyz_xz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_0_yz[i] = -2.0 * g_z_xz_0_yz[i] * a_exp + 4.0 * g_yyz_xz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_0_zz[i] = -2.0 * g_z_xz_0_zz[i] * a_exp + 4.0 * g_yyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_yy_0_0_0_z_yy_0_xx, g_yy_0_0_0_z_yy_0_xy, g_yy_0_0_0_z_yy_0_xz, g_yy_0_0_0_z_yy_0_yy, g_yy_0_0_0_z_yy_0_yz, g_yy_0_0_0_z_yy_0_zz, g_yyz_yy_0_xx, g_yyz_yy_0_xy, g_yyz_yy_0_xz, g_yyz_yy_0_yy, g_yyz_yy_0_yz, g_yyz_yy_0_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_yy_0_xx[i] = -2.0 * g_z_yy_0_xx[i] * a_exp + 4.0 * g_yyz_yy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_0_xy[i] = -2.0 * g_z_yy_0_xy[i] * a_exp + 4.0 * g_yyz_yy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_0_xz[i] = -2.0 * g_z_yy_0_xz[i] * a_exp + 4.0 * g_yyz_yy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_0_yy[i] = -2.0 * g_z_yy_0_yy[i] * a_exp + 4.0 * g_yyz_yy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_0_yz[i] = -2.0 * g_z_yy_0_yz[i] * a_exp + 4.0 * g_yyz_yy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_0_zz[i] = -2.0 * g_z_yy_0_zz[i] * a_exp + 4.0 * g_yyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_yy_0_0_0_z_yz_0_xx, g_yy_0_0_0_z_yz_0_xy, g_yy_0_0_0_z_yz_0_xz, g_yy_0_0_0_z_yz_0_yy, g_yy_0_0_0_z_yz_0_yz, g_yy_0_0_0_z_yz_0_zz, g_yyz_yz_0_xx, g_yyz_yz_0_xy, g_yyz_yz_0_xz, g_yyz_yz_0_yy, g_yyz_yz_0_yz, g_yyz_yz_0_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_yz_0_xx[i] = -2.0 * g_z_yz_0_xx[i] * a_exp + 4.0 * g_yyz_yz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_0_xy[i] = -2.0 * g_z_yz_0_xy[i] * a_exp + 4.0 * g_yyz_yz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_0_xz[i] = -2.0 * g_z_yz_0_xz[i] * a_exp + 4.0 * g_yyz_yz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_0_yy[i] = -2.0 * g_z_yz_0_yy[i] * a_exp + 4.0 * g_yyz_yz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_0_yz[i] = -2.0 * g_z_yz_0_yz[i] * a_exp + 4.0 * g_yyz_yz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_0_zz[i] = -2.0 * g_z_yz_0_zz[i] * a_exp + 4.0 * g_yyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_yy_0_0_0_z_zz_0_xx, g_yy_0_0_0_z_zz_0_xy, g_yy_0_0_0_z_zz_0_xz, g_yy_0_0_0_z_zz_0_yy, g_yy_0_0_0_z_zz_0_yz, g_yy_0_0_0_z_zz_0_zz, g_yyz_zz_0_xx, g_yyz_zz_0_xy, g_yyz_zz_0_xz, g_yyz_zz_0_yy, g_yyz_zz_0_yz, g_yyz_zz_0_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_zz_0_xx[i] = -2.0 * g_z_zz_0_xx[i] * a_exp + 4.0 * g_yyz_zz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_0_xy[i] = -2.0 * g_z_zz_0_xy[i] * a_exp + 4.0 * g_yyz_zz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_0_xz[i] = -2.0 * g_z_zz_0_xz[i] * a_exp + 4.0 * g_yyz_zz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_0_yy[i] = -2.0 * g_z_zz_0_yy[i] * a_exp + 4.0 * g_yyz_zz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_0_yz[i] = -2.0 * g_z_zz_0_yz[i] * a_exp + 4.0 * g_yyz_zz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_0_zz[i] = -2.0 * g_z_zz_0_zz[i] * a_exp + 4.0 * g_yyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_xyz_xx_0_xx, g_xyz_xx_0_xy, g_xyz_xx_0_xz, g_xyz_xx_0_yy, g_xyz_xx_0_yz, g_xyz_xx_0_zz, g_yz_0_0_0_x_xx_0_xx, g_yz_0_0_0_x_xx_0_xy, g_yz_0_0_0_x_xx_0_xz, g_yz_0_0_0_x_xx_0_yy, g_yz_0_0_0_x_xx_0_yz, g_yz_0_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xx_0_xx[i] = 4.0 * g_xyz_xx_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_0_xy[i] = 4.0 * g_xyz_xx_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_0_xz[i] = 4.0 * g_xyz_xx_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_0_yy[i] = 4.0 * g_xyz_xx_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_0_yz[i] = 4.0 * g_xyz_xx_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_0_zz[i] = 4.0 * g_xyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_xyz_xy_0_xx, g_xyz_xy_0_xy, g_xyz_xy_0_xz, g_xyz_xy_0_yy, g_xyz_xy_0_yz, g_xyz_xy_0_zz, g_yz_0_0_0_x_xy_0_xx, g_yz_0_0_0_x_xy_0_xy, g_yz_0_0_0_x_xy_0_xz, g_yz_0_0_0_x_xy_0_yy, g_yz_0_0_0_x_xy_0_yz, g_yz_0_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xy_0_xx[i] = 4.0 * g_xyz_xy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_0_xy[i] = 4.0 * g_xyz_xy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_0_xz[i] = 4.0 * g_xyz_xy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_0_yy[i] = 4.0 * g_xyz_xy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_0_yz[i] = 4.0 * g_xyz_xy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_0_zz[i] = 4.0 * g_xyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_xyz_xz_0_xx, g_xyz_xz_0_xy, g_xyz_xz_0_xz, g_xyz_xz_0_yy, g_xyz_xz_0_yz, g_xyz_xz_0_zz, g_yz_0_0_0_x_xz_0_xx, g_yz_0_0_0_x_xz_0_xy, g_yz_0_0_0_x_xz_0_xz, g_yz_0_0_0_x_xz_0_yy, g_yz_0_0_0_x_xz_0_yz, g_yz_0_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xz_0_xx[i] = 4.0 * g_xyz_xz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_0_xy[i] = 4.0 * g_xyz_xz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_0_xz[i] = 4.0 * g_xyz_xz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_0_yy[i] = 4.0 * g_xyz_xz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_0_yz[i] = 4.0 * g_xyz_xz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_0_zz[i] = 4.0 * g_xyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_xyz_yy_0_xx, g_xyz_yy_0_xy, g_xyz_yy_0_xz, g_xyz_yy_0_yy, g_xyz_yy_0_yz, g_xyz_yy_0_zz, g_yz_0_0_0_x_yy_0_xx, g_yz_0_0_0_x_yy_0_xy, g_yz_0_0_0_x_yy_0_xz, g_yz_0_0_0_x_yy_0_yy, g_yz_0_0_0_x_yy_0_yz, g_yz_0_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_yy_0_xx[i] = 4.0 * g_xyz_yy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_0_xy[i] = 4.0 * g_xyz_yy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_0_xz[i] = 4.0 * g_xyz_yy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_0_yy[i] = 4.0 * g_xyz_yy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_0_yz[i] = 4.0 * g_xyz_yy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_0_zz[i] = 4.0 * g_xyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_xyz_yz_0_xx, g_xyz_yz_0_xy, g_xyz_yz_0_xz, g_xyz_yz_0_yy, g_xyz_yz_0_yz, g_xyz_yz_0_zz, g_yz_0_0_0_x_yz_0_xx, g_yz_0_0_0_x_yz_0_xy, g_yz_0_0_0_x_yz_0_xz, g_yz_0_0_0_x_yz_0_yy, g_yz_0_0_0_x_yz_0_yz, g_yz_0_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_yz_0_xx[i] = 4.0 * g_xyz_yz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_0_xy[i] = 4.0 * g_xyz_yz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_0_xz[i] = 4.0 * g_xyz_yz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_0_yy[i] = 4.0 * g_xyz_yz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_0_yz[i] = 4.0 * g_xyz_yz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_0_zz[i] = 4.0 * g_xyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_xyz_zz_0_xx, g_xyz_zz_0_xy, g_xyz_zz_0_xz, g_xyz_zz_0_yy, g_xyz_zz_0_yz, g_xyz_zz_0_zz, g_yz_0_0_0_x_zz_0_xx, g_yz_0_0_0_x_zz_0_xy, g_yz_0_0_0_x_zz_0_xz, g_yz_0_0_0_x_zz_0_yy, g_yz_0_0_0_x_zz_0_yz, g_yz_0_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_zz_0_xx[i] = 4.0 * g_xyz_zz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_0_xy[i] = 4.0 * g_xyz_zz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_0_xz[i] = 4.0 * g_xyz_zz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_0_yy[i] = 4.0 * g_xyz_zz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_0_yz[i] = 4.0 * g_xyz_zz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_0_zz[i] = 4.0 * g_xyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_yyz_xx_0_xx, g_yyz_xx_0_xy, g_yyz_xx_0_xz, g_yyz_xx_0_yy, g_yyz_xx_0_yz, g_yyz_xx_0_zz, g_yz_0_0_0_y_xx_0_xx, g_yz_0_0_0_y_xx_0_xy, g_yz_0_0_0_y_xx_0_xz, g_yz_0_0_0_y_xx_0_yy, g_yz_0_0_0_y_xx_0_yz, g_yz_0_0_0_y_xx_0_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xx_0_xx[i] = -2.0 * g_z_xx_0_xx[i] * a_exp + 4.0 * g_yyz_xx_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_0_xy[i] = -2.0 * g_z_xx_0_xy[i] * a_exp + 4.0 * g_yyz_xx_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_0_xz[i] = -2.0 * g_z_xx_0_xz[i] * a_exp + 4.0 * g_yyz_xx_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_0_yy[i] = -2.0 * g_z_xx_0_yy[i] * a_exp + 4.0 * g_yyz_xx_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_0_yz[i] = -2.0 * g_z_xx_0_yz[i] * a_exp + 4.0 * g_yyz_xx_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_0_zz[i] = -2.0 * g_z_xx_0_zz[i] * a_exp + 4.0 * g_yyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_yyz_xy_0_xx, g_yyz_xy_0_xy, g_yyz_xy_0_xz, g_yyz_xy_0_yy, g_yyz_xy_0_yz, g_yyz_xy_0_zz, g_yz_0_0_0_y_xy_0_xx, g_yz_0_0_0_y_xy_0_xy, g_yz_0_0_0_y_xy_0_xz, g_yz_0_0_0_y_xy_0_yy, g_yz_0_0_0_y_xy_0_yz, g_yz_0_0_0_y_xy_0_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xy_0_xx[i] = -2.0 * g_z_xy_0_xx[i] * a_exp + 4.0 * g_yyz_xy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_0_xy[i] = -2.0 * g_z_xy_0_xy[i] * a_exp + 4.0 * g_yyz_xy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_0_xz[i] = -2.0 * g_z_xy_0_xz[i] * a_exp + 4.0 * g_yyz_xy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_0_yy[i] = -2.0 * g_z_xy_0_yy[i] * a_exp + 4.0 * g_yyz_xy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_0_yz[i] = -2.0 * g_z_xy_0_yz[i] * a_exp + 4.0 * g_yyz_xy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_0_zz[i] = -2.0 * g_z_xy_0_zz[i] * a_exp + 4.0 * g_yyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_yyz_xz_0_xx, g_yyz_xz_0_xy, g_yyz_xz_0_xz, g_yyz_xz_0_yy, g_yyz_xz_0_yz, g_yyz_xz_0_zz, g_yz_0_0_0_y_xz_0_xx, g_yz_0_0_0_y_xz_0_xy, g_yz_0_0_0_y_xz_0_xz, g_yz_0_0_0_y_xz_0_yy, g_yz_0_0_0_y_xz_0_yz, g_yz_0_0_0_y_xz_0_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xz_0_xx[i] = -2.0 * g_z_xz_0_xx[i] * a_exp + 4.0 * g_yyz_xz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_0_xy[i] = -2.0 * g_z_xz_0_xy[i] * a_exp + 4.0 * g_yyz_xz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_0_xz[i] = -2.0 * g_z_xz_0_xz[i] * a_exp + 4.0 * g_yyz_xz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_0_yy[i] = -2.0 * g_z_xz_0_yy[i] * a_exp + 4.0 * g_yyz_xz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_0_yz[i] = -2.0 * g_z_xz_0_yz[i] * a_exp + 4.0 * g_yyz_xz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_0_zz[i] = -2.0 * g_z_xz_0_zz[i] * a_exp + 4.0 * g_yyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_yyz_yy_0_xx, g_yyz_yy_0_xy, g_yyz_yy_0_xz, g_yyz_yy_0_yy, g_yyz_yy_0_yz, g_yyz_yy_0_zz, g_yz_0_0_0_y_yy_0_xx, g_yz_0_0_0_y_yy_0_xy, g_yz_0_0_0_y_yy_0_xz, g_yz_0_0_0_y_yy_0_yy, g_yz_0_0_0_y_yy_0_yz, g_yz_0_0_0_y_yy_0_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_yy_0_xx[i] = -2.0 * g_z_yy_0_xx[i] * a_exp + 4.0 * g_yyz_yy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_0_xy[i] = -2.0 * g_z_yy_0_xy[i] * a_exp + 4.0 * g_yyz_yy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_0_xz[i] = -2.0 * g_z_yy_0_xz[i] * a_exp + 4.0 * g_yyz_yy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_0_yy[i] = -2.0 * g_z_yy_0_yy[i] * a_exp + 4.0 * g_yyz_yy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_0_yz[i] = -2.0 * g_z_yy_0_yz[i] * a_exp + 4.0 * g_yyz_yy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_0_zz[i] = -2.0 * g_z_yy_0_zz[i] * a_exp + 4.0 * g_yyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_yyz_yz_0_xx, g_yyz_yz_0_xy, g_yyz_yz_0_xz, g_yyz_yz_0_yy, g_yyz_yz_0_yz, g_yyz_yz_0_zz, g_yz_0_0_0_y_yz_0_xx, g_yz_0_0_0_y_yz_0_xy, g_yz_0_0_0_y_yz_0_xz, g_yz_0_0_0_y_yz_0_yy, g_yz_0_0_0_y_yz_0_yz, g_yz_0_0_0_y_yz_0_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_yz_0_xx[i] = -2.0 * g_z_yz_0_xx[i] * a_exp + 4.0 * g_yyz_yz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_0_xy[i] = -2.0 * g_z_yz_0_xy[i] * a_exp + 4.0 * g_yyz_yz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_0_xz[i] = -2.0 * g_z_yz_0_xz[i] * a_exp + 4.0 * g_yyz_yz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_0_yy[i] = -2.0 * g_z_yz_0_yy[i] * a_exp + 4.0 * g_yyz_yz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_0_yz[i] = -2.0 * g_z_yz_0_yz[i] * a_exp + 4.0 * g_yyz_yz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_0_zz[i] = -2.0 * g_z_yz_0_zz[i] * a_exp + 4.0 * g_yyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_yyz_zz_0_xx, g_yyz_zz_0_xy, g_yyz_zz_0_xz, g_yyz_zz_0_yy, g_yyz_zz_0_yz, g_yyz_zz_0_zz, g_yz_0_0_0_y_zz_0_xx, g_yz_0_0_0_y_zz_0_xy, g_yz_0_0_0_y_zz_0_xz, g_yz_0_0_0_y_zz_0_yy, g_yz_0_0_0_y_zz_0_yz, g_yz_0_0_0_y_zz_0_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_zz_0_xx[i] = -2.0 * g_z_zz_0_xx[i] * a_exp + 4.0 * g_yyz_zz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_0_xy[i] = -2.0 * g_z_zz_0_xy[i] * a_exp + 4.0 * g_yyz_zz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_0_xz[i] = -2.0 * g_z_zz_0_xz[i] * a_exp + 4.0 * g_yyz_zz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_0_yy[i] = -2.0 * g_z_zz_0_yy[i] * a_exp + 4.0 * g_yyz_zz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_0_yz[i] = -2.0 * g_z_zz_0_yz[i] * a_exp + 4.0 * g_yyz_zz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_0_zz[i] = -2.0 * g_z_zz_0_zz[i] * a_exp + 4.0 * g_yyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz, g_yz_0_0_0_z_xx_0_xx, g_yz_0_0_0_z_xx_0_xy, g_yz_0_0_0_z_xx_0_xz, g_yz_0_0_0_z_xx_0_yy, g_yz_0_0_0_z_xx_0_yz, g_yz_0_0_0_z_xx_0_zz, g_yzz_xx_0_xx, g_yzz_xx_0_xy, g_yzz_xx_0_xz, g_yzz_xx_0_yy, g_yzz_xx_0_yz, g_yzz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xx_0_xx[i] = -2.0 * g_y_xx_0_xx[i] * a_exp + 4.0 * g_yzz_xx_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_0_xy[i] = -2.0 * g_y_xx_0_xy[i] * a_exp + 4.0 * g_yzz_xx_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_0_xz[i] = -2.0 * g_y_xx_0_xz[i] * a_exp + 4.0 * g_yzz_xx_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_0_yy[i] = -2.0 * g_y_xx_0_yy[i] * a_exp + 4.0 * g_yzz_xx_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_0_yz[i] = -2.0 * g_y_xx_0_yz[i] * a_exp + 4.0 * g_yzz_xx_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_0_zz[i] = -2.0 * g_y_xx_0_zz[i] * a_exp + 4.0 * g_yzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz, g_yz_0_0_0_z_xy_0_xx, g_yz_0_0_0_z_xy_0_xy, g_yz_0_0_0_z_xy_0_xz, g_yz_0_0_0_z_xy_0_yy, g_yz_0_0_0_z_xy_0_yz, g_yz_0_0_0_z_xy_0_zz, g_yzz_xy_0_xx, g_yzz_xy_0_xy, g_yzz_xy_0_xz, g_yzz_xy_0_yy, g_yzz_xy_0_yz, g_yzz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xy_0_xx[i] = -2.0 * g_y_xy_0_xx[i] * a_exp + 4.0 * g_yzz_xy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_0_xy[i] = -2.0 * g_y_xy_0_xy[i] * a_exp + 4.0 * g_yzz_xy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_0_xz[i] = -2.0 * g_y_xy_0_xz[i] * a_exp + 4.0 * g_yzz_xy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_0_yy[i] = -2.0 * g_y_xy_0_yy[i] * a_exp + 4.0 * g_yzz_xy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_0_yz[i] = -2.0 * g_y_xy_0_yz[i] * a_exp + 4.0 * g_yzz_xy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_0_zz[i] = -2.0 * g_y_xy_0_zz[i] * a_exp + 4.0 * g_yzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz, g_yz_0_0_0_z_xz_0_xx, g_yz_0_0_0_z_xz_0_xy, g_yz_0_0_0_z_xz_0_xz, g_yz_0_0_0_z_xz_0_yy, g_yz_0_0_0_z_xz_0_yz, g_yz_0_0_0_z_xz_0_zz, g_yzz_xz_0_xx, g_yzz_xz_0_xy, g_yzz_xz_0_xz, g_yzz_xz_0_yy, g_yzz_xz_0_yz, g_yzz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xz_0_xx[i] = -2.0 * g_y_xz_0_xx[i] * a_exp + 4.0 * g_yzz_xz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_0_xy[i] = -2.0 * g_y_xz_0_xy[i] * a_exp + 4.0 * g_yzz_xz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_0_xz[i] = -2.0 * g_y_xz_0_xz[i] * a_exp + 4.0 * g_yzz_xz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_0_yy[i] = -2.0 * g_y_xz_0_yy[i] * a_exp + 4.0 * g_yzz_xz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_0_yz[i] = -2.0 * g_y_xz_0_yz[i] * a_exp + 4.0 * g_yzz_xz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_0_zz[i] = -2.0 * g_y_xz_0_zz[i] * a_exp + 4.0 * g_yzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz, g_yz_0_0_0_z_yy_0_xx, g_yz_0_0_0_z_yy_0_xy, g_yz_0_0_0_z_yy_0_xz, g_yz_0_0_0_z_yy_0_yy, g_yz_0_0_0_z_yy_0_yz, g_yz_0_0_0_z_yy_0_zz, g_yzz_yy_0_xx, g_yzz_yy_0_xy, g_yzz_yy_0_xz, g_yzz_yy_0_yy, g_yzz_yy_0_yz, g_yzz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_yy_0_xx[i] = -2.0 * g_y_yy_0_xx[i] * a_exp + 4.0 * g_yzz_yy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_0_xy[i] = -2.0 * g_y_yy_0_xy[i] * a_exp + 4.0 * g_yzz_yy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_0_xz[i] = -2.0 * g_y_yy_0_xz[i] * a_exp + 4.0 * g_yzz_yy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_0_yy[i] = -2.0 * g_y_yy_0_yy[i] * a_exp + 4.0 * g_yzz_yy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_0_yz[i] = -2.0 * g_y_yy_0_yz[i] * a_exp + 4.0 * g_yzz_yy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_0_zz[i] = -2.0 * g_y_yy_0_zz[i] * a_exp + 4.0 * g_yzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz, g_yz_0_0_0_z_yz_0_xx, g_yz_0_0_0_z_yz_0_xy, g_yz_0_0_0_z_yz_0_xz, g_yz_0_0_0_z_yz_0_yy, g_yz_0_0_0_z_yz_0_yz, g_yz_0_0_0_z_yz_0_zz, g_yzz_yz_0_xx, g_yzz_yz_0_xy, g_yzz_yz_0_xz, g_yzz_yz_0_yy, g_yzz_yz_0_yz, g_yzz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_yz_0_xx[i] = -2.0 * g_y_yz_0_xx[i] * a_exp + 4.0 * g_yzz_yz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_0_xy[i] = -2.0 * g_y_yz_0_xy[i] * a_exp + 4.0 * g_yzz_yz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_0_xz[i] = -2.0 * g_y_yz_0_xz[i] * a_exp + 4.0 * g_yzz_yz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_0_yy[i] = -2.0 * g_y_yz_0_yy[i] * a_exp + 4.0 * g_yzz_yz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_0_yz[i] = -2.0 * g_y_yz_0_yz[i] * a_exp + 4.0 * g_yzz_yz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_0_zz[i] = -2.0 * g_y_yz_0_zz[i] * a_exp + 4.0 * g_yzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz, g_yz_0_0_0_z_zz_0_xx, g_yz_0_0_0_z_zz_0_xy, g_yz_0_0_0_z_zz_0_xz, g_yz_0_0_0_z_zz_0_yy, g_yz_0_0_0_z_zz_0_yz, g_yz_0_0_0_z_zz_0_zz, g_yzz_zz_0_xx, g_yzz_zz_0_xy, g_yzz_zz_0_xz, g_yzz_zz_0_yy, g_yzz_zz_0_yz, g_yzz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_zz_0_xx[i] = -2.0 * g_y_zz_0_xx[i] * a_exp + 4.0 * g_yzz_zz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_0_xy[i] = -2.0 * g_y_zz_0_xy[i] * a_exp + 4.0 * g_yzz_zz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_0_xz[i] = -2.0 * g_y_zz_0_xz[i] * a_exp + 4.0 * g_yzz_zz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_0_yy[i] = -2.0 * g_y_zz_0_yy[i] * a_exp + 4.0 * g_yzz_zz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_0_yz[i] = -2.0 * g_y_zz_0_yz[i] * a_exp + 4.0 * g_yzz_zz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_0_zz[i] = -2.0 * g_y_zz_0_zz[i] * a_exp + 4.0 * g_yzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_xzz_xx_0_xx, g_xzz_xx_0_xy, g_xzz_xx_0_xz, g_xzz_xx_0_yy, g_xzz_xx_0_yz, g_xzz_xx_0_zz, g_zz_0_0_0_x_xx_0_xx, g_zz_0_0_0_x_xx_0_xy, g_zz_0_0_0_x_xx_0_xz, g_zz_0_0_0_x_xx_0_yy, g_zz_0_0_0_x_xx_0_yz, g_zz_0_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xx_0_xx[i] = -2.0 * g_x_xx_0_xx[i] * a_exp + 4.0 * g_xzz_xx_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_0_xy[i] = -2.0 * g_x_xx_0_xy[i] * a_exp + 4.0 * g_xzz_xx_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_0_xz[i] = -2.0 * g_x_xx_0_xz[i] * a_exp + 4.0 * g_xzz_xx_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_0_yy[i] = -2.0 * g_x_xx_0_yy[i] * a_exp + 4.0 * g_xzz_xx_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_0_yz[i] = -2.0 * g_x_xx_0_yz[i] * a_exp + 4.0 * g_xzz_xx_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_0_zz[i] = -2.0 * g_x_xx_0_zz[i] * a_exp + 4.0 * g_xzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_xzz_xy_0_xx, g_xzz_xy_0_xy, g_xzz_xy_0_xz, g_xzz_xy_0_yy, g_xzz_xy_0_yz, g_xzz_xy_0_zz, g_zz_0_0_0_x_xy_0_xx, g_zz_0_0_0_x_xy_0_xy, g_zz_0_0_0_x_xy_0_xz, g_zz_0_0_0_x_xy_0_yy, g_zz_0_0_0_x_xy_0_yz, g_zz_0_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xy_0_xx[i] = -2.0 * g_x_xy_0_xx[i] * a_exp + 4.0 * g_xzz_xy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_0_xy[i] = -2.0 * g_x_xy_0_xy[i] * a_exp + 4.0 * g_xzz_xy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_0_xz[i] = -2.0 * g_x_xy_0_xz[i] * a_exp + 4.0 * g_xzz_xy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_0_yy[i] = -2.0 * g_x_xy_0_yy[i] * a_exp + 4.0 * g_xzz_xy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_0_yz[i] = -2.0 * g_x_xy_0_yz[i] * a_exp + 4.0 * g_xzz_xy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_0_zz[i] = -2.0 * g_x_xy_0_zz[i] * a_exp + 4.0 * g_xzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_xzz_xz_0_xx, g_xzz_xz_0_xy, g_xzz_xz_0_xz, g_xzz_xz_0_yy, g_xzz_xz_0_yz, g_xzz_xz_0_zz, g_zz_0_0_0_x_xz_0_xx, g_zz_0_0_0_x_xz_0_xy, g_zz_0_0_0_x_xz_0_xz, g_zz_0_0_0_x_xz_0_yy, g_zz_0_0_0_x_xz_0_yz, g_zz_0_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xz_0_xx[i] = -2.0 * g_x_xz_0_xx[i] * a_exp + 4.0 * g_xzz_xz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_0_xy[i] = -2.0 * g_x_xz_0_xy[i] * a_exp + 4.0 * g_xzz_xz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_0_xz[i] = -2.0 * g_x_xz_0_xz[i] * a_exp + 4.0 * g_xzz_xz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_0_yy[i] = -2.0 * g_x_xz_0_yy[i] * a_exp + 4.0 * g_xzz_xz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_0_yz[i] = -2.0 * g_x_xz_0_yz[i] * a_exp + 4.0 * g_xzz_xz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_0_zz[i] = -2.0 * g_x_xz_0_zz[i] * a_exp + 4.0 * g_xzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_xzz_yy_0_xx, g_xzz_yy_0_xy, g_xzz_yy_0_xz, g_xzz_yy_0_yy, g_xzz_yy_0_yz, g_xzz_yy_0_zz, g_zz_0_0_0_x_yy_0_xx, g_zz_0_0_0_x_yy_0_xy, g_zz_0_0_0_x_yy_0_xz, g_zz_0_0_0_x_yy_0_yy, g_zz_0_0_0_x_yy_0_yz, g_zz_0_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_yy_0_xx[i] = -2.0 * g_x_yy_0_xx[i] * a_exp + 4.0 * g_xzz_yy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_0_xy[i] = -2.0 * g_x_yy_0_xy[i] * a_exp + 4.0 * g_xzz_yy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_0_xz[i] = -2.0 * g_x_yy_0_xz[i] * a_exp + 4.0 * g_xzz_yy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_0_yy[i] = -2.0 * g_x_yy_0_yy[i] * a_exp + 4.0 * g_xzz_yy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_0_yz[i] = -2.0 * g_x_yy_0_yz[i] * a_exp + 4.0 * g_xzz_yy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_0_zz[i] = -2.0 * g_x_yy_0_zz[i] * a_exp + 4.0 * g_xzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_xzz_yz_0_xx, g_xzz_yz_0_xy, g_xzz_yz_0_xz, g_xzz_yz_0_yy, g_xzz_yz_0_yz, g_xzz_yz_0_zz, g_zz_0_0_0_x_yz_0_xx, g_zz_0_0_0_x_yz_0_xy, g_zz_0_0_0_x_yz_0_xz, g_zz_0_0_0_x_yz_0_yy, g_zz_0_0_0_x_yz_0_yz, g_zz_0_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_yz_0_xx[i] = -2.0 * g_x_yz_0_xx[i] * a_exp + 4.0 * g_xzz_yz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_0_xy[i] = -2.0 * g_x_yz_0_xy[i] * a_exp + 4.0 * g_xzz_yz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_0_xz[i] = -2.0 * g_x_yz_0_xz[i] * a_exp + 4.0 * g_xzz_yz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_0_yy[i] = -2.0 * g_x_yz_0_yy[i] * a_exp + 4.0 * g_xzz_yz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_0_yz[i] = -2.0 * g_x_yz_0_yz[i] * a_exp + 4.0 * g_xzz_yz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_0_zz[i] = -2.0 * g_x_yz_0_zz[i] * a_exp + 4.0 * g_xzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_xzz_zz_0_xx, g_xzz_zz_0_xy, g_xzz_zz_0_xz, g_xzz_zz_0_yy, g_xzz_zz_0_yz, g_xzz_zz_0_zz, g_zz_0_0_0_x_zz_0_xx, g_zz_0_0_0_x_zz_0_xy, g_zz_0_0_0_x_zz_0_xz, g_zz_0_0_0_x_zz_0_yy, g_zz_0_0_0_x_zz_0_yz, g_zz_0_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_zz_0_xx[i] = -2.0 * g_x_zz_0_xx[i] * a_exp + 4.0 * g_xzz_zz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_0_xy[i] = -2.0 * g_x_zz_0_xy[i] * a_exp + 4.0 * g_xzz_zz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_0_xz[i] = -2.0 * g_x_zz_0_xz[i] * a_exp + 4.0 * g_xzz_zz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_0_yy[i] = -2.0 * g_x_zz_0_yy[i] * a_exp + 4.0 * g_xzz_zz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_0_yz[i] = -2.0 * g_x_zz_0_yz[i] * a_exp + 4.0 * g_xzz_zz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_0_zz[i] = -2.0 * g_x_zz_0_zz[i] * a_exp + 4.0 * g_xzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz, g_yzz_xx_0_xx, g_yzz_xx_0_xy, g_yzz_xx_0_xz, g_yzz_xx_0_yy, g_yzz_xx_0_yz, g_yzz_xx_0_zz, g_zz_0_0_0_y_xx_0_xx, g_zz_0_0_0_y_xx_0_xy, g_zz_0_0_0_y_xx_0_xz, g_zz_0_0_0_y_xx_0_yy, g_zz_0_0_0_y_xx_0_yz, g_zz_0_0_0_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xx_0_xx[i] = -2.0 * g_y_xx_0_xx[i] * a_exp + 4.0 * g_yzz_xx_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_0_xy[i] = -2.0 * g_y_xx_0_xy[i] * a_exp + 4.0 * g_yzz_xx_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_0_xz[i] = -2.0 * g_y_xx_0_xz[i] * a_exp + 4.0 * g_yzz_xx_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_0_yy[i] = -2.0 * g_y_xx_0_yy[i] * a_exp + 4.0 * g_yzz_xx_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_0_yz[i] = -2.0 * g_y_xx_0_yz[i] * a_exp + 4.0 * g_yzz_xx_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_0_zz[i] = -2.0 * g_y_xx_0_zz[i] * a_exp + 4.0 * g_yzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz, g_yzz_xy_0_xx, g_yzz_xy_0_xy, g_yzz_xy_0_xz, g_yzz_xy_0_yy, g_yzz_xy_0_yz, g_yzz_xy_0_zz, g_zz_0_0_0_y_xy_0_xx, g_zz_0_0_0_y_xy_0_xy, g_zz_0_0_0_y_xy_0_xz, g_zz_0_0_0_y_xy_0_yy, g_zz_0_0_0_y_xy_0_yz, g_zz_0_0_0_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xy_0_xx[i] = -2.0 * g_y_xy_0_xx[i] * a_exp + 4.0 * g_yzz_xy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_0_xy[i] = -2.0 * g_y_xy_0_xy[i] * a_exp + 4.0 * g_yzz_xy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_0_xz[i] = -2.0 * g_y_xy_0_xz[i] * a_exp + 4.0 * g_yzz_xy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_0_yy[i] = -2.0 * g_y_xy_0_yy[i] * a_exp + 4.0 * g_yzz_xy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_0_yz[i] = -2.0 * g_y_xy_0_yz[i] * a_exp + 4.0 * g_yzz_xy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_0_zz[i] = -2.0 * g_y_xy_0_zz[i] * a_exp + 4.0 * g_yzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz, g_yzz_xz_0_xx, g_yzz_xz_0_xy, g_yzz_xz_0_xz, g_yzz_xz_0_yy, g_yzz_xz_0_yz, g_yzz_xz_0_zz, g_zz_0_0_0_y_xz_0_xx, g_zz_0_0_0_y_xz_0_xy, g_zz_0_0_0_y_xz_0_xz, g_zz_0_0_0_y_xz_0_yy, g_zz_0_0_0_y_xz_0_yz, g_zz_0_0_0_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xz_0_xx[i] = -2.0 * g_y_xz_0_xx[i] * a_exp + 4.0 * g_yzz_xz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_0_xy[i] = -2.0 * g_y_xz_0_xy[i] * a_exp + 4.0 * g_yzz_xz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_0_xz[i] = -2.0 * g_y_xz_0_xz[i] * a_exp + 4.0 * g_yzz_xz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_0_yy[i] = -2.0 * g_y_xz_0_yy[i] * a_exp + 4.0 * g_yzz_xz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_0_yz[i] = -2.0 * g_y_xz_0_yz[i] * a_exp + 4.0 * g_yzz_xz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_0_zz[i] = -2.0 * g_y_xz_0_zz[i] * a_exp + 4.0 * g_yzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz, g_yzz_yy_0_xx, g_yzz_yy_0_xy, g_yzz_yy_0_xz, g_yzz_yy_0_yy, g_yzz_yy_0_yz, g_yzz_yy_0_zz, g_zz_0_0_0_y_yy_0_xx, g_zz_0_0_0_y_yy_0_xy, g_zz_0_0_0_y_yy_0_xz, g_zz_0_0_0_y_yy_0_yy, g_zz_0_0_0_y_yy_0_yz, g_zz_0_0_0_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_yy_0_xx[i] = -2.0 * g_y_yy_0_xx[i] * a_exp + 4.0 * g_yzz_yy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_0_xy[i] = -2.0 * g_y_yy_0_xy[i] * a_exp + 4.0 * g_yzz_yy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_0_xz[i] = -2.0 * g_y_yy_0_xz[i] * a_exp + 4.0 * g_yzz_yy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_0_yy[i] = -2.0 * g_y_yy_0_yy[i] * a_exp + 4.0 * g_yzz_yy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_0_yz[i] = -2.0 * g_y_yy_0_yz[i] * a_exp + 4.0 * g_yzz_yy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_0_zz[i] = -2.0 * g_y_yy_0_zz[i] * a_exp + 4.0 * g_yzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz, g_yzz_yz_0_xx, g_yzz_yz_0_xy, g_yzz_yz_0_xz, g_yzz_yz_0_yy, g_yzz_yz_0_yz, g_yzz_yz_0_zz, g_zz_0_0_0_y_yz_0_xx, g_zz_0_0_0_y_yz_0_xy, g_zz_0_0_0_y_yz_0_xz, g_zz_0_0_0_y_yz_0_yy, g_zz_0_0_0_y_yz_0_yz, g_zz_0_0_0_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_yz_0_xx[i] = -2.0 * g_y_yz_0_xx[i] * a_exp + 4.0 * g_yzz_yz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_0_xy[i] = -2.0 * g_y_yz_0_xy[i] * a_exp + 4.0 * g_yzz_yz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_0_xz[i] = -2.0 * g_y_yz_0_xz[i] * a_exp + 4.0 * g_yzz_yz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_0_yy[i] = -2.0 * g_y_yz_0_yy[i] * a_exp + 4.0 * g_yzz_yz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_0_yz[i] = -2.0 * g_y_yz_0_yz[i] * a_exp + 4.0 * g_yzz_yz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_0_zz[i] = -2.0 * g_y_yz_0_zz[i] * a_exp + 4.0 * g_yzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz, g_yzz_zz_0_xx, g_yzz_zz_0_xy, g_yzz_zz_0_xz, g_yzz_zz_0_yy, g_yzz_zz_0_yz, g_yzz_zz_0_zz, g_zz_0_0_0_y_zz_0_xx, g_zz_0_0_0_y_zz_0_xy, g_zz_0_0_0_y_zz_0_xz, g_zz_0_0_0_y_zz_0_yy, g_zz_0_0_0_y_zz_0_yz, g_zz_0_0_0_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_zz_0_xx[i] = -2.0 * g_y_zz_0_xx[i] * a_exp + 4.0 * g_yzz_zz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_0_xy[i] = -2.0 * g_y_zz_0_xy[i] * a_exp + 4.0 * g_yzz_zz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_0_xz[i] = -2.0 * g_y_zz_0_xz[i] * a_exp + 4.0 * g_yzz_zz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_0_yy[i] = -2.0 * g_y_zz_0_yy[i] * a_exp + 4.0 * g_yzz_zz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_0_yz[i] = -2.0 * g_y_zz_0_yz[i] * a_exp + 4.0 * g_yzz_zz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_0_zz[i] = -2.0 * g_y_zz_0_zz[i] * a_exp + 4.0 * g_yzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz, g_zz_0_0_0_z_xx_0_xx, g_zz_0_0_0_z_xx_0_xy, g_zz_0_0_0_z_xx_0_xz, g_zz_0_0_0_z_xx_0_yy, g_zz_0_0_0_z_xx_0_yz, g_zz_0_0_0_z_xx_0_zz, g_zzz_xx_0_xx, g_zzz_xx_0_xy, g_zzz_xx_0_xz, g_zzz_xx_0_yy, g_zzz_xx_0_yz, g_zzz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xx_0_xx[i] = -6.0 * g_z_xx_0_xx[i] * a_exp + 4.0 * g_zzz_xx_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_0_xy[i] = -6.0 * g_z_xx_0_xy[i] * a_exp + 4.0 * g_zzz_xx_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_0_xz[i] = -6.0 * g_z_xx_0_xz[i] * a_exp + 4.0 * g_zzz_xx_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_0_yy[i] = -6.0 * g_z_xx_0_yy[i] * a_exp + 4.0 * g_zzz_xx_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_0_yz[i] = -6.0 * g_z_xx_0_yz[i] * a_exp + 4.0 * g_zzz_xx_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_0_zz[i] = -6.0 * g_z_xx_0_zz[i] * a_exp + 4.0 * g_zzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz, g_zz_0_0_0_z_xy_0_xx, g_zz_0_0_0_z_xy_0_xy, g_zz_0_0_0_z_xy_0_xz, g_zz_0_0_0_z_xy_0_yy, g_zz_0_0_0_z_xy_0_yz, g_zz_0_0_0_z_xy_0_zz, g_zzz_xy_0_xx, g_zzz_xy_0_xy, g_zzz_xy_0_xz, g_zzz_xy_0_yy, g_zzz_xy_0_yz, g_zzz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xy_0_xx[i] = -6.0 * g_z_xy_0_xx[i] * a_exp + 4.0 * g_zzz_xy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_0_xy[i] = -6.0 * g_z_xy_0_xy[i] * a_exp + 4.0 * g_zzz_xy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_0_xz[i] = -6.0 * g_z_xy_0_xz[i] * a_exp + 4.0 * g_zzz_xy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_0_yy[i] = -6.0 * g_z_xy_0_yy[i] * a_exp + 4.0 * g_zzz_xy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_0_yz[i] = -6.0 * g_z_xy_0_yz[i] * a_exp + 4.0 * g_zzz_xy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_0_zz[i] = -6.0 * g_z_xy_0_zz[i] * a_exp + 4.0 * g_zzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz, g_zz_0_0_0_z_xz_0_xx, g_zz_0_0_0_z_xz_0_xy, g_zz_0_0_0_z_xz_0_xz, g_zz_0_0_0_z_xz_0_yy, g_zz_0_0_0_z_xz_0_yz, g_zz_0_0_0_z_xz_0_zz, g_zzz_xz_0_xx, g_zzz_xz_0_xy, g_zzz_xz_0_xz, g_zzz_xz_0_yy, g_zzz_xz_0_yz, g_zzz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xz_0_xx[i] = -6.0 * g_z_xz_0_xx[i] * a_exp + 4.0 * g_zzz_xz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_0_xy[i] = -6.0 * g_z_xz_0_xy[i] * a_exp + 4.0 * g_zzz_xz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_0_xz[i] = -6.0 * g_z_xz_0_xz[i] * a_exp + 4.0 * g_zzz_xz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_0_yy[i] = -6.0 * g_z_xz_0_yy[i] * a_exp + 4.0 * g_zzz_xz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_0_yz[i] = -6.0 * g_z_xz_0_yz[i] * a_exp + 4.0 * g_zzz_xz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_0_zz[i] = -6.0 * g_z_xz_0_zz[i] * a_exp + 4.0 * g_zzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz, g_zz_0_0_0_z_yy_0_xx, g_zz_0_0_0_z_yy_0_xy, g_zz_0_0_0_z_yy_0_xz, g_zz_0_0_0_z_yy_0_yy, g_zz_0_0_0_z_yy_0_yz, g_zz_0_0_0_z_yy_0_zz, g_zzz_yy_0_xx, g_zzz_yy_0_xy, g_zzz_yy_0_xz, g_zzz_yy_0_yy, g_zzz_yy_0_yz, g_zzz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_yy_0_xx[i] = -6.0 * g_z_yy_0_xx[i] * a_exp + 4.0 * g_zzz_yy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_0_xy[i] = -6.0 * g_z_yy_0_xy[i] * a_exp + 4.0 * g_zzz_yy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_0_xz[i] = -6.0 * g_z_yy_0_xz[i] * a_exp + 4.0 * g_zzz_yy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_0_yy[i] = -6.0 * g_z_yy_0_yy[i] * a_exp + 4.0 * g_zzz_yy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_0_yz[i] = -6.0 * g_z_yy_0_yz[i] * a_exp + 4.0 * g_zzz_yy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_0_zz[i] = -6.0 * g_z_yy_0_zz[i] * a_exp + 4.0 * g_zzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz, g_zz_0_0_0_z_yz_0_xx, g_zz_0_0_0_z_yz_0_xy, g_zz_0_0_0_z_yz_0_xz, g_zz_0_0_0_z_yz_0_yy, g_zz_0_0_0_z_yz_0_yz, g_zz_0_0_0_z_yz_0_zz, g_zzz_yz_0_xx, g_zzz_yz_0_xy, g_zzz_yz_0_xz, g_zzz_yz_0_yy, g_zzz_yz_0_yz, g_zzz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_yz_0_xx[i] = -6.0 * g_z_yz_0_xx[i] * a_exp + 4.0 * g_zzz_yz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_0_xy[i] = -6.0 * g_z_yz_0_xy[i] * a_exp + 4.0 * g_zzz_yz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_0_xz[i] = -6.0 * g_z_yz_0_xz[i] * a_exp + 4.0 * g_zzz_yz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_0_yy[i] = -6.0 * g_z_yz_0_yy[i] * a_exp + 4.0 * g_zzz_yz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_0_yz[i] = -6.0 * g_z_yz_0_yz[i] * a_exp + 4.0 * g_zzz_yz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_0_zz[i] = -6.0 * g_z_yz_0_zz[i] * a_exp + 4.0 * g_zzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz, g_zz_0_0_0_z_zz_0_xx, g_zz_0_0_0_z_zz_0_xy, g_zz_0_0_0_z_zz_0_xz, g_zz_0_0_0_z_zz_0_yy, g_zz_0_0_0_z_zz_0_yz, g_zz_0_0_0_z_zz_0_zz, g_zzz_zz_0_xx, g_zzz_zz_0_xy, g_zzz_zz_0_xz, g_zzz_zz_0_yy, g_zzz_zz_0_yz, g_zzz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_zz_0_xx[i] = -6.0 * g_z_zz_0_xx[i] * a_exp + 4.0 * g_zzz_zz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_0_xy[i] = -6.0 * g_z_zz_0_xy[i] * a_exp + 4.0 * g_zzz_zz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_0_xz[i] = -6.0 * g_z_zz_0_xz[i] * a_exp + 4.0 * g_zzz_zz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_0_yy[i] = -6.0 * g_z_zz_0_yy[i] * a_exp + 4.0 * g_zzz_zz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_0_yz[i] = -6.0 * g_z_zz_0_yz[i] * a_exp + 4.0 * g_zzz_zz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_0_zz[i] = -6.0 * g_z_zz_0_zz[i] * a_exp + 4.0 * g_zzz_zz_0_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

