#include "GeomDeriv1000OfScalarForDDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_ddsd_0(CSimdArray<double>& buffer_1000_ddsd,
                     const CSimdArray<double>& buffer_pdsd,
                     const CSimdArray<double>& buffer_fdsd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_ddsd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_ddsd

    auto g_x_0_0_0_xx_xx_0_xx = buffer_1000_ddsd[0];

    auto g_x_0_0_0_xx_xx_0_xy = buffer_1000_ddsd[1];

    auto g_x_0_0_0_xx_xx_0_xz = buffer_1000_ddsd[2];

    auto g_x_0_0_0_xx_xx_0_yy = buffer_1000_ddsd[3];

    auto g_x_0_0_0_xx_xx_0_yz = buffer_1000_ddsd[4];

    auto g_x_0_0_0_xx_xx_0_zz = buffer_1000_ddsd[5];

    auto g_x_0_0_0_xx_xy_0_xx = buffer_1000_ddsd[6];

    auto g_x_0_0_0_xx_xy_0_xy = buffer_1000_ddsd[7];

    auto g_x_0_0_0_xx_xy_0_xz = buffer_1000_ddsd[8];

    auto g_x_0_0_0_xx_xy_0_yy = buffer_1000_ddsd[9];

    auto g_x_0_0_0_xx_xy_0_yz = buffer_1000_ddsd[10];

    auto g_x_0_0_0_xx_xy_0_zz = buffer_1000_ddsd[11];

    auto g_x_0_0_0_xx_xz_0_xx = buffer_1000_ddsd[12];

    auto g_x_0_0_0_xx_xz_0_xy = buffer_1000_ddsd[13];

    auto g_x_0_0_0_xx_xz_0_xz = buffer_1000_ddsd[14];

    auto g_x_0_0_0_xx_xz_0_yy = buffer_1000_ddsd[15];

    auto g_x_0_0_0_xx_xz_0_yz = buffer_1000_ddsd[16];

    auto g_x_0_0_0_xx_xz_0_zz = buffer_1000_ddsd[17];

    auto g_x_0_0_0_xx_yy_0_xx = buffer_1000_ddsd[18];

    auto g_x_0_0_0_xx_yy_0_xy = buffer_1000_ddsd[19];

    auto g_x_0_0_0_xx_yy_0_xz = buffer_1000_ddsd[20];

    auto g_x_0_0_0_xx_yy_0_yy = buffer_1000_ddsd[21];

    auto g_x_0_0_0_xx_yy_0_yz = buffer_1000_ddsd[22];

    auto g_x_0_0_0_xx_yy_0_zz = buffer_1000_ddsd[23];

    auto g_x_0_0_0_xx_yz_0_xx = buffer_1000_ddsd[24];

    auto g_x_0_0_0_xx_yz_0_xy = buffer_1000_ddsd[25];

    auto g_x_0_0_0_xx_yz_0_xz = buffer_1000_ddsd[26];

    auto g_x_0_0_0_xx_yz_0_yy = buffer_1000_ddsd[27];

    auto g_x_0_0_0_xx_yz_0_yz = buffer_1000_ddsd[28];

    auto g_x_0_0_0_xx_yz_0_zz = buffer_1000_ddsd[29];

    auto g_x_0_0_0_xx_zz_0_xx = buffer_1000_ddsd[30];

    auto g_x_0_0_0_xx_zz_0_xy = buffer_1000_ddsd[31];

    auto g_x_0_0_0_xx_zz_0_xz = buffer_1000_ddsd[32];

    auto g_x_0_0_0_xx_zz_0_yy = buffer_1000_ddsd[33];

    auto g_x_0_0_0_xx_zz_0_yz = buffer_1000_ddsd[34];

    auto g_x_0_0_0_xx_zz_0_zz = buffer_1000_ddsd[35];

    auto g_x_0_0_0_xy_xx_0_xx = buffer_1000_ddsd[36];

    auto g_x_0_0_0_xy_xx_0_xy = buffer_1000_ddsd[37];

    auto g_x_0_0_0_xy_xx_0_xz = buffer_1000_ddsd[38];

    auto g_x_0_0_0_xy_xx_0_yy = buffer_1000_ddsd[39];

    auto g_x_0_0_0_xy_xx_0_yz = buffer_1000_ddsd[40];

    auto g_x_0_0_0_xy_xx_0_zz = buffer_1000_ddsd[41];

    auto g_x_0_0_0_xy_xy_0_xx = buffer_1000_ddsd[42];

    auto g_x_0_0_0_xy_xy_0_xy = buffer_1000_ddsd[43];

    auto g_x_0_0_0_xy_xy_0_xz = buffer_1000_ddsd[44];

    auto g_x_0_0_0_xy_xy_0_yy = buffer_1000_ddsd[45];

    auto g_x_0_0_0_xy_xy_0_yz = buffer_1000_ddsd[46];

    auto g_x_0_0_0_xy_xy_0_zz = buffer_1000_ddsd[47];

    auto g_x_0_0_0_xy_xz_0_xx = buffer_1000_ddsd[48];

    auto g_x_0_0_0_xy_xz_0_xy = buffer_1000_ddsd[49];

    auto g_x_0_0_0_xy_xz_0_xz = buffer_1000_ddsd[50];

    auto g_x_0_0_0_xy_xz_0_yy = buffer_1000_ddsd[51];

    auto g_x_0_0_0_xy_xz_0_yz = buffer_1000_ddsd[52];

    auto g_x_0_0_0_xy_xz_0_zz = buffer_1000_ddsd[53];

    auto g_x_0_0_0_xy_yy_0_xx = buffer_1000_ddsd[54];

    auto g_x_0_0_0_xy_yy_0_xy = buffer_1000_ddsd[55];

    auto g_x_0_0_0_xy_yy_0_xz = buffer_1000_ddsd[56];

    auto g_x_0_0_0_xy_yy_0_yy = buffer_1000_ddsd[57];

    auto g_x_0_0_0_xy_yy_0_yz = buffer_1000_ddsd[58];

    auto g_x_0_0_0_xy_yy_0_zz = buffer_1000_ddsd[59];

    auto g_x_0_0_0_xy_yz_0_xx = buffer_1000_ddsd[60];

    auto g_x_0_0_0_xy_yz_0_xy = buffer_1000_ddsd[61];

    auto g_x_0_0_0_xy_yz_0_xz = buffer_1000_ddsd[62];

    auto g_x_0_0_0_xy_yz_0_yy = buffer_1000_ddsd[63];

    auto g_x_0_0_0_xy_yz_0_yz = buffer_1000_ddsd[64];

    auto g_x_0_0_0_xy_yz_0_zz = buffer_1000_ddsd[65];

    auto g_x_0_0_0_xy_zz_0_xx = buffer_1000_ddsd[66];

    auto g_x_0_0_0_xy_zz_0_xy = buffer_1000_ddsd[67];

    auto g_x_0_0_0_xy_zz_0_xz = buffer_1000_ddsd[68];

    auto g_x_0_0_0_xy_zz_0_yy = buffer_1000_ddsd[69];

    auto g_x_0_0_0_xy_zz_0_yz = buffer_1000_ddsd[70];

    auto g_x_0_0_0_xy_zz_0_zz = buffer_1000_ddsd[71];

    auto g_x_0_0_0_xz_xx_0_xx = buffer_1000_ddsd[72];

    auto g_x_0_0_0_xz_xx_0_xy = buffer_1000_ddsd[73];

    auto g_x_0_0_0_xz_xx_0_xz = buffer_1000_ddsd[74];

    auto g_x_0_0_0_xz_xx_0_yy = buffer_1000_ddsd[75];

    auto g_x_0_0_0_xz_xx_0_yz = buffer_1000_ddsd[76];

    auto g_x_0_0_0_xz_xx_0_zz = buffer_1000_ddsd[77];

    auto g_x_0_0_0_xz_xy_0_xx = buffer_1000_ddsd[78];

    auto g_x_0_0_0_xz_xy_0_xy = buffer_1000_ddsd[79];

    auto g_x_0_0_0_xz_xy_0_xz = buffer_1000_ddsd[80];

    auto g_x_0_0_0_xz_xy_0_yy = buffer_1000_ddsd[81];

    auto g_x_0_0_0_xz_xy_0_yz = buffer_1000_ddsd[82];

    auto g_x_0_0_0_xz_xy_0_zz = buffer_1000_ddsd[83];

    auto g_x_0_0_0_xz_xz_0_xx = buffer_1000_ddsd[84];

    auto g_x_0_0_0_xz_xz_0_xy = buffer_1000_ddsd[85];

    auto g_x_0_0_0_xz_xz_0_xz = buffer_1000_ddsd[86];

    auto g_x_0_0_0_xz_xz_0_yy = buffer_1000_ddsd[87];

    auto g_x_0_0_0_xz_xz_0_yz = buffer_1000_ddsd[88];

    auto g_x_0_0_0_xz_xz_0_zz = buffer_1000_ddsd[89];

    auto g_x_0_0_0_xz_yy_0_xx = buffer_1000_ddsd[90];

    auto g_x_0_0_0_xz_yy_0_xy = buffer_1000_ddsd[91];

    auto g_x_0_0_0_xz_yy_0_xz = buffer_1000_ddsd[92];

    auto g_x_0_0_0_xz_yy_0_yy = buffer_1000_ddsd[93];

    auto g_x_0_0_0_xz_yy_0_yz = buffer_1000_ddsd[94];

    auto g_x_0_0_0_xz_yy_0_zz = buffer_1000_ddsd[95];

    auto g_x_0_0_0_xz_yz_0_xx = buffer_1000_ddsd[96];

    auto g_x_0_0_0_xz_yz_0_xy = buffer_1000_ddsd[97];

    auto g_x_0_0_0_xz_yz_0_xz = buffer_1000_ddsd[98];

    auto g_x_0_0_0_xz_yz_0_yy = buffer_1000_ddsd[99];

    auto g_x_0_0_0_xz_yz_0_yz = buffer_1000_ddsd[100];

    auto g_x_0_0_0_xz_yz_0_zz = buffer_1000_ddsd[101];

    auto g_x_0_0_0_xz_zz_0_xx = buffer_1000_ddsd[102];

    auto g_x_0_0_0_xz_zz_0_xy = buffer_1000_ddsd[103];

    auto g_x_0_0_0_xz_zz_0_xz = buffer_1000_ddsd[104];

    auto g_x_0_0_0_xz_zz_0_yy = buffer_1000_ddsd[105];

    auto g_x_0_0_0_xz_zz_0_yz = buffer_1000_ddsd[106];

    auto g_x_0_0_0_xz_zz_0_zz = buffer_1000_ddsd[107];

    auto g_x_0_0_0_yy_xx_0_xx = buffer_1000_ddsd[108];

    auto g_x_0_0_0_yy_xx_0_xy = buffer_1000_ddsd[109];

    auto g_x_0_0_0_yy_xx_0_xz = buffer_1000_ddsd[110];

    auto g_x_0_0_0_yy_xx_0_yy = buffer_1000_ddsd[111];

    auto g_x_0_0_0_yy_xx_0_yz = buffer_1000_ddsd[112];

    auto g_x_0_0_0_yy_xx_0_zz = buffer_1000_ddsd[113];

    auto g_x_0_0_0_yy_xy_0_xx = buffer_1000_ddsd[114];

    auto g_x_0_0_0_yy_xy_0_xy = buffer_1000_ddsd[115];

    auto g_x_0_0_0_yy_xy_0_xz = buffer_1000_ddsd[116];

    auto g_x_0_0_0_yy_xy_0_yy = buffer_1000_ddsd[117];

    auto g_x_0_0_0_yy_xy_0_yz = buffer_1000_ddsd[118];

    auto g_x_0_0_0_yy_xy_0_zz = buffer_1000_ddsd[119];

    auto g_x_0_0_0_yy_xz_0_xx = buffer_1000_ddsd[120];

    auto g_x_0_0_0_yy_xz_0_xy = buffer_1000_ddsd[121];

    auto g_x_0_0_0_yy_xz_0_xz = buffer_1000_ddsd[122];

    auto g_x_0_0_0_yy_xz_0_yy = buffer_1000_ddsd[123];

    auto g_x_0_0_0_yy_xz_0_yz = buffer_1000_ddsd[124];

    auto g_x_0_0_0_yy_xz_0_zz = buffer_1000_ddsd[125];

    auto g_x_0_0_0_yy_yy_0_xx = buffer_1000_ddsd[126];

    auto g_x_0_0_0_yy_yy_0_xy = buffer_1000_ddsd[127];

    auto g_x_0_0_0_yy_yy_0_xz = buffer_1000_ddsd[128];

    auto g_x_0_0_0_yy_yy_0_yy = buffer_1000_ddsd[129];

    auto g_x_0_0_0_yy_yy_0_yz = buffer_1000_ddsd[130];

    auto g_x_0_0_0_yy_yy_0_zz = buffer_1000_ddsd[131];

    auto g_x_0_0_0_yy_yz_0_xx = buffer_1000_ddsd[132];

    auto g_x_0_0_0_yy_yz_0_xy = buffer_1000_ddsd[133];

    auto g_x_0_0_0_yy_yz_0_xz = buffer_1000_ddsd[134];

    auto g_x_0_0_0_yy_yz_0_yy = buffer_1000_ddsd[135];

    auto g_x_0_0_0_yy_yz_0_yz = buffer_1000_ddsd[136];

    auto g_x_0_0_0_yy_yz_0_zz = buffer_1000_ddsd[137];

    auto g_x_0_0_0_yy_zz_0_xx = buffer_1000_ddsd[138];

    auto g_x_0_0_0_yy_zz_0_xy = buffer_1000_ddsd[139];

    auto g_x_0_0_0_yy_zz_0_xz = buffer_1000_ddsd[140];

    auto g_x_0_0_0_yy_zz_0_yy = buffer_1000_ddsd[141];

    auto g_x_0_0_0_yy_zz_0_yz = buffer_1000_ddsd[142];

    auto g_x_0_0_0_yy_zz_0_zz = buffer_1000_ddsd[143];

    auto g_x_0_0_0_yz_xx_0_xx = buffer_1000_ddsd[144];

    auto g_x_0_0_0_yz_xx_0_xy = buffer_1000_ddsd[145];

    auto g_x_0_0_0_yz_xx_0_xz = buffer_1000_ddsd[146];

    auto g_x_0_0_0_yz_xx_0_yy = buffer_1000_ddsd[147];

    auto g_x_0_0_0_yz_xx_0_yz = buffer_1000_ddsd[148];

    auto g_x_0_0_0_yz_xx_0_zz = buffer_1000_ddsd[149];

    auto g_x_0_0_0_yz_xy_0_xx = buffer_1000_ddsd[150];

    auto g_x_0_0_0_yz_xy_0_xy = buffer_1000_ddsd[151];

    auto g_x_0_0_0_yz_xy_0_xz = buffer_1000_ddsd[152];

    auto g_x_0_0_0_yz_xy_0_yy = buffer_1000_ddsd[153];

    auto g_x_0_0_0_yz_xy_0_yz = buffer_1000_ddsd[154];

    auto g_x_0_0_0_yz_xy_0_zz = buffer_1000_ddsd[155];

    auto g_x_0_0_0_yz_xz_0_xx = buffer_1000_ddsd[156];

    auto g_x_0_0_0_yz_xz_0_xy = buffer_1000_ddsd[157];

    auto g_x_0_0_0_yz_xz_0_xz = buffer_1000_ddsd[158];

    auto g_x_0_0_0_yz_xz_0_yy = buffer_1000_ddsd[159];

    auto g_x_0_0_0_yz_xz_0_yz = buffer_1000_ddsd[160];

    auto g_x_0_0_0_yz_xz_0_zz = buffer_1000_ddsd[161];

    auto g_x_0_0_0_yz_yy_0_xx = buffer_1000_ddsd[162];

    auto g_x_0_0_0_yz_yy_0_xy = buffer_1000_ddsd[163];

    auto g_x_0_0_0_yz_yy_0_xz = buffer_1000_ddsd[164];

    auto g_x_0_0_0_yz_yy_0_yy = buffer_1000_ddsd[165];

    auto g_x_0_0_0_yz_yy_0_yz = buffer_1000_ddsd[166];

    auto g_x_0_0_0_yz_yy_0_zz = buffer_1000_ddsd[167];

    auto g_x_0_0_0_yz_yz_0_xx = buffer_1000_ddsd[168];

    auto g_x_0_0_0_yz_yz_0_xy = buffer_1000_ddsd[169];

    auto g_x_0_0_0_yz_yz_0_xz = buffer_1000_ddsd[170];

    auto g_x_0_0_0_yz_yz_0_yy = buffer_1000_ddsd[171];

    auto g_x_0_0_0_yz_yz_0_yz = buffer_1000_ddsd[172];

    auto g_x_0_0_0_yz_yz_0_zz = buffer_1000_ddsd[173];

    auto g_x_0_0_0_yz_zz_0_xx = buffer_1000_ddsd[174];

    auto g_x_0_0_0_yz_zz_0_xy = buffer_1000_ddsd[175];

    auto g_x_0_0_0_yz_zz_0_xz = buffer_1000_ddsd[176];

    auto g_x_0_0_0_yz_zz_0_yy = buffer_1000_ddsd[177];

    auto g_x_0_0_0_yz_zz_0_yz = buffer_1000_ddsd[178];

    auto g_x_0_0_0_yz_zz_0_zz = buffer_1000_ddsd[179];

    auto g_x_0_0_0_zz_xx_0_xx = buffer_1000_ddsd[180];

    auto g_x_0_0_0_zz_xx_0_xy = buffer_1000_ddsd[181];

    auto g_x_0_0_0_zz_xx_0_xz = buffer_1000_ddsd[182];

    auto g_x_0_0_0_zz_xx_0_yy = buffer_1000_ddsd[183];

    auto g_x_0_0_0_zz_xx_0_yz = buffer_1000_ddsd[184];

    auto g_x_0_0_0_zz_xx_0_zz = buffer_1000_ddsd[185];

    auto g_x_0_0_0_zz_xy_0_xx = buffer_1000_ddsd[186];

    auto g_x_0_0_0_zz_xy_0_xy = buffer_1000_ddsd[187];

    auto g_x_0_0_0_zz_xy_0_xz = buffer_1000_ddsd[188];

    auto g_x_0_0_0_zz_xy_0_yy = buffer_1000_ddsd[189];

    auto g_x_0_0_0_zz_xy_0_yz = buffer_1000_ddsd[190];

    auto g_x_0_0_0_zz_xy_0_zz = buffer_1000_ddsd[191];

    auto g_x_0_0_0_zz_xz_0_xx = buffer_1000_ddsd[192];

    auto g_x_0_0_0_zz_xz_0_xy = buffer_1000_ddsd[193];

    auto g_x_0_0_0_zz_xz_0_xz = buffer_1000_ddsd[194];

    auto g_x_0_0_0_zz_xz_0_yy = buffer_1000_ddsd[195];

    auto g_x_0_0_0_zz_xz_0_yz = buffer_1000_ddsd[196];

    auto g_x_0_0_0_zz_xz_0_zz = buffer_1000_ddsd[197];

    auto g_x_0_0_0_zz_yy_0_xx = buffer_1000_ddsd[198];

    auto g_x_0_0_0_zz_yy_0_xy = buffer_1000_ddsd[199];

    auto g_x_0_0_0_zz_yy_0_xz = buffer_1000_ddsd[200];

    auto g_x_0_0_0_zz_yy_0_yy = buffer_1000_ddsd[201];

    auto g_x_0_0_0_zz_yy_0_yz = buffer_1000_ddsd[202];

    auto g_x_0_0_0_zz_yy_0_zz = buffer_1000_ddsd[203];

    auto g_x_0_0_0_zz_yz_0_xx = buffer_1000_ddsd[204];

    auto g_x_0_0_0_zz_yz_0_xy = buffer_1000_ddsd[205];

    auto g_x_0_0_0_zz_yz_0_xz = buffer_1000_ddsd[206];

    auto g_x_0_0_0_zz_yz_0_yy = buffer_1000_ddsd[207];

    auto g_x_0_0_0_zz_yz_0_yz = buffer_1000_ddsd[208];

    auto g_x_0_0_0_zz_yz_0_zz = buffer_1000_ddsd[209];

    auto g_x_0_0_0_zz_zz_0_xx = buffer_1000_ddsd[210];

    auto g_x_0_0_0_zz_zz_0_xy = buffer_1000_ddsd[211];

    auto g_x_0_0_0_zz_zz_0_xz = buffer_1000_ddsd[212];

    auto g_x_0_0_0_zz_zz_0_yy = buffer_1000_ddsd[213];

    auto g_x_0_0_0_zz_zz_0_yz = buffer_1000_ddsd[214];

    auto g_x_0_0_0_zz_zz_0_zz = buffer_1000_ddsd[215];

    auto g_y_0_0_0_xx_xx_0_xx = buffer_1000_ddsd[216];

    auto g_y_0_0_0_xx_xx_0_xy = buffer_1000_ddsd[217];

    auto g_y_0_0_0_xx_xx_0_xz = buffer_1000_ddsd[218];

    auto g_y_0_0_0_xx_xx_0_yy = buffer_1000_ddsd[219];

    auto g_y_0_0_0_xx_xx_0_yz = buffer_1000_ddsd[220];

    auto g_y_0_0_0_xx_xx_0_zz = buffer_1000_ddsd[221];

    auto g_y_0_0_0_xx_xy_0_xx = buffer_1000_ddsd[222];

    auto g_y_0_0_0_xx_xy_0_xy = buffer_1000_ddsd[223];

    auto g_y_0_0_0_xx_xy_0_xz = buffer_1000_ddsd[224];

    auto g_y_0_0_0_xx_xy_0_yy = buffer_1000_ddsd[225];

    auto g_y_0_0_0_xx_xy_0_yz = buffer_1000_ddsd[226];

    auto g_y_0_0_0_xx_xy_0_zz = buffer_1000_ddsd[227];

    auto g_y_0_0_0_xx_xz_0_xx = buffer_1000_ddsd[228];

    auto g_y_0_0_0_xx_xz_0_xy = buffer_1000_ddsd[229];

    auto g_y_0_0_0_xx_xz_0_xz = buffer_1000_ddsd[230];

    auto g_y_0_0_0_xx_xz_0_yy = buffer_1000_ddsd[231];

    auto g_y_0_0_0_xx_xz_0_yz = buffer_1000_ddsd[232];

    auto g_y_0_0_0_xx_xz_0_zz = buffer_1000_ddsd[233];

    auto g_y_0_0_0_xx_yy_0_xx = buffer_1000_ddsd[234];

    auto g_y_0_0_0_xx_yy_0_xy = buffer_1000_ddsd[235];

    auto g_y_0_0_0_xx_yy_0_xz = buffer_1000_ddsd[236];

    auto g_y_0_0_0_xx_yy_0_yy = buffer_1000_ddsd[237];

    auto g_y_0_0_0_xx_yy_0_yz = buffer_1000_ddsd[238];

    auto g_y_0_0_0_xx_yy_0_zz = buffer_1000_ddsd[239];

    auto g_y_0_0_0_xx_yz_0_xx = buffer_1000_ddsd[240];

    auto g_y_0_0_0_xx_yz_0_xy = buffer_1000_ddsd[241];

    auto g_y_0_0_0_xx_yz_0_xz = buffer_1000_ddsd[242];

    auto g_y_0_0_0_xx_yz_0_yy = buffer_1000_ddsd[243];

    auto g_y_0_0_0_xx_yz_0_yz = buffer_1000_ddsd[244];

    auto g_y_0_0_0_xx_yz_0_zz = buffer_1000_ddsd[245];

    auto g_y_0_0_0_xx_zz_0_xx = buffer_1000_ddsd[246];

    auto g_y_0_0_0_xx_zz_0_xy = buffer_1000_ddsd[247];

    auto g_y_0_0_0_xx_zz_0_xz = buffer_1000_ddsd[248];

    auto g_y_0_0_0_xx_zz_0_yy = buffer_1000_ddsd[249];

    auto g_y_0_0_0_xx_zz_0_yz = buffer_1000_ddsd[250];

    auto g_y_0_0_0_xx_zz_0_zz = buffer_1000_ddsd[251];

    auto g_y_0_0_0_xy_xx_0_xx = buffer_1000_ddsd[252];

    auto g_y_0_0_0_xy_xx_0_xy = buffer_1000_ddsd[253];

    auto g_y_0_0_0_xy_xx_0_xz = buffer_1000_ddsd[254];

    auto g_y_0_0_0_xy_xx_0_yy = buffer_1000_ddsd[255];

    auto g_y_0_0_0_xy_xx_0_yz = buffer_1000_ddsd[256];

    auto g_y_0_0_0_xy_xx_0_zz = buffer_1000_ddsd[257];

    auto g_y_0_0_0_xy_xy_0_xx = buffer_1000_ddsd[258];

    auto g_y_0_0_0_xy_xy_0_xy = buffer_1000_ddsd[259];

    auto g_y_0_0_0_xy_xy_0_xz = buffer_1000_ddsd[260];

    auto g_y_0_0_0_xy_xy_0_yy = buffer_1000_ddsd[261];

    auto g_y_0_0_0_xy_xy_0_yz = buffer_1000_ddsd[262];

    auto g_y_0_0_0_xy_xy_0_zz = buffer_1000_ddsd[263];

    auto g_y_0_0_0_xy_xz_0_xx = buffer_1000_ddsd[264];

    auto g_y_0_0_0_xy_xz_0_xy = buffer_1000_ddsd[265];

    auto g_y_0_0_0_xy_xz_0_xz = buffer_1000_ddsd[266];

    auto g_y_0_0_0_xy_xz_0_yy = buffer_1000_ddsd[267];

    auto g_y_0_0_0_xy_xz_0_yz = buffer_1000_ddsd[268];

    auto g_y_0_0_0_xy_xz_0_zz = buffer_1000_ddsd[269];

    auto g_y_0_0_0_xy_yy_0_xx = buffer_1000_ddsd[270];

    auto g_y_0_0_0_xy_yy_0_xy = buffer_1000_ddsd[271];

    auto g_y_0_0_0_xy_yy_0_xz = buffer_1000_ddsd[272];

    auto g_y_0_0_0_xy_yy_0_yy = buffer_1000_ddsd[273];

    auto g_y_0_0_0_xy_yy_0_yz = buffer_1000_ddsd[274];

    auto g_y_0_0_0_xy_yy_0_zz = buffer_1000_ddsd[275];

    auto g_y_0_0_0_xy_yz_0_xx = buffer_1000_ddsd[276];

    auto g_y_0_0_0_xy_yz_0_xy = buffer_1000_ddsd[277];

    auto g_y_0_0_0_xy_yz_0_xz = buffer_1000_ddsd[278];

    auto g_y_0_0_0_xy_yz_0_yy = buffer_1000_ddsd[279];

    auto g_y_0_0_0_xy_yz_0_yz = buffer_1000_ddsd[280];

    auto g_y_0_0_0_xy_yz_0_zz = buffer_1000_ddsd[281];

    auto g_y_0_0_0_xy_zz_0_xx = buffer_1000_ddsd[282];

    auto g_y_0_0_0_xy_zz_0_xy = buffer_1000_ddsd[283];

    auto g_y_0_0_0_xy_zz_0_xz = buffer_1000_ddsd[284];

    auto g_y_0_0_0_xy_zz_0_yy = buffer_1000_ddsd[285];

    auto g_y_0_0_0_xy_zz_0_yz = buffer_1000_ddsd[286];

    auto g_y_0_0_0_xy_zz_0_zz = buffer_1000_ddsd[287];

    auto g_y_0_0_0_xz_xx_0_xx = buffer_1000_ddsd[288];

    auto g_y_0_0_0_xz_xx_0_xy = buffer_1000_ddsd[289];

    auto g_y_0_0_0_xz_xx_0_xz = buffer_1000_ddsd[290];

    auto g_y_0_0_0_xz_xx_0_yy = buffer_1000_ddsd[291];

    auto g_y_0_0_0_xz_xx_0_yz = buffer_1000_ddsd[292];

    auto g_y_0_0_0_xz_xx_0_zz = buffer_1000_ddsd[293];

    auto g_y_0_0_0_xz_xy_0_xx = buffer_1000_ddsd[294];

    auto g_y_0_0_0_xz_xy_0_xy = buffer_1000_ddsd[295];

    auto g_y_0_0_0_xz_xy_0_xz = buffer_1000_ddsd[296];

    auto g_y_0_0_0_xz_xy_0_yy = buffer_1000_ddsd[297];

    auto g_y_0_0_0_xz_xy_0_yz = buffer_1000_ddsd[298];

    auto g_y_0_0_0_xz_xy_0_zz = buffer_1000_ddsd[299];

    auto g_y_0_0_0_xz_xz_0_xx = buffer_1000_ddsd[300];

    auto g_y_0_0_0_xz_xz_0_xy = buffer_1000_ddsd[301];

    auto g_y_0_0_0_xz_xz_0_xz = buffer_1000_ddsd[302];

    auto g_y_0_0_0_xz_xz_0_yy = buffer_1000_ddsd[303];

    auto g_y_0_0_0_xz_xz_0_yz = buffer_1000_ddsd[304];

    auto g_y_0_0_0_xz_xz_0_zz = buffer_1000_ddsd[305];

    auto g_y_0_0_0_xz_yy_0_xx = buffer_1000_ddsd[306];

    auto g_y_0_0_0_xz_yy_0_xy = buffer_1000_ddsd[307];

    auto g_y_0_0_0_xz_yy_0_xz = buffer_1000_ddsd[308];

    auto g_y_0_0_0_xz_yy_0_yy = buffer_1000_ddsd[309];

    auto g_y_0_0_0_xz_yy_0_yz = buffer_1000_ddsd[310];

    auto g_y_0_0_0_xz_yy_0_zz = buffer_1000_ddsd[311];

    auto g_y_0_0_0_xz_yz_0_xx = buffer_1000_ddsd[312];

    auto g_y_0_0_0_xz_yz_0_xy = buffer_1000_ddsd[313];

    auto g_y_0_0_0_xz_yz_0_xz = buffer_1000_ddsd[314];

    auto g_y_0_0_0_xz_yz_0_yy = buffer_1000_ddsd[315];

    auto g_y_0_0_0_xz_yz_0_yz = buffer_1000_ddsd[316];

    auto g_y_0_0_0_xz_yz_0_zz = buffer_1000_ddsd[317];

    auto g_y_0_0_0_xz_zz_0_xx = buffer_1000_ddsd[318];

    auto g_y_0_0_0_xz_zz_0_xy = buffer_1000_ddsd[319];

    auto g_y_0_0_0_xz_zz_0_xz = buffer_1000_ddsd[320];

    auto g_y_0_0_0_xz_zz_0_yy = buffer_1000_ddsd[321];

    auto g_y_0_0_0_xz_zz_0_yz = buffer_1000_ddsd[322];

    auto g_y_0_0_0_xz_zz_0_zz = buffer_1000_ddsd[323];

    auto g_y_0_0_0_yy_xx_0_xx = buffer_1000_ddsd[324];

    auto g_y_0_0_0_yy_xx_0_xy = buffer_1000_ddsd[325];

    auto g_y_0_0_0_yy_xx_0_xz = buffer_1000_ddsd[326];

    auto g_y_0_0_0_yy_xx_0_yy = buffer_1000_ddsd[327];

    auto g_y_0_0_0_yy_xx_0_yz = buffer_1000_ddsd[328];

    auto g_y_0_0_0_yy_xx_0_zz = buffer_1000_ddsd[329];

    auto g_y_0_0_0_yy_xy_0_xx = buffer_1000_ddsd[330];

    auto g_y_0_0_0_yy_xy_0_xy = buffer_1000_ddsd[331];

    auto g_y_0_0_0_yy_xy_0_xz = buffer_1000_ddsd[332];

    auto g_y_0_0_0_yy_xy_0_yy = buffer_1000_ddsd[333];

    auto g_y_0_0_0_yy_xy_0_yz = buffer_1000_ddsd[334];

    auto g_y_0_0_0_yy_xy_0_zz = buffer_1000_ddsd[335];

    auto g_y_0_0_0_yy_xz_0_xx = buffer_1000_ddsd[336];

    auto g_y_0_0_0_yy_xz_0_xy = buffer_1000_ddsd[337];

    auto g_y_0_0_0_yy_xz_0_xz = buffer_1000_ddsd[338];

    auto g_y_0_0_0_yy_xz_0_yy = buffer_1000_ddsd[339];

    auto g_y_0_0_0_yy_xz_0_yz = buffer_1000_ddsd[340];

    auto g_y_0_0_0_yy_xz_0_zz = buffer_1000_ddsd[341];

    auto g_y_0_0_0_yy_yy_0_xx = buffer_1000_ddsd[342];

    auto g_y_0_0_0_yy_yy_0_xy = buffer_1000_ddsd[343];

    auto g_y_0_0_0_yy_yy_0_xz = buffer_1000_ddsd[344];

    auto g_y_0_0_0_yy_yy_0_yy = buffer_1000_ddsd[345];

    auto g_y_0_0_0_yy_yy_0_yz = buffer_1000_ddsd[346];

    auto g_y_0_0_0_yy_yy_0_zz = buffer_1000_ddsd[347];

    auto g_y_0_0_0_yy_yz_0_xx = buffer_1000_ddsd[348];

    auto g_y_0_0_0_yy_yz_0_xy = buffer_1000_ddsd[349];

    auto g_y_0_0_0_yy_yz_0_xz = buffer_1000_ddsd[350];

    auto g_y_0_0_0_yy_yz_0_yy = buffer_1000_ddsd[351];

    auto g_y_0_0_0_yy_yz_0_yz = buffer_1000_ddsd[352];

    auto g_y_0_0_0_yy_yz_0_zz = buffer_1000_ddsd[353];

    auto g_y_0_0_0_yy_zz_0_xx = buffer_1000_ddsd[354];

    auto g_y_0_0_0_yy_zz_0_xy = buffer_1000_ddsd[355];

    auto g_y_0_0_0_yy_zz_0_xz = buffer_1000_ddsd[356];

    auto g_y_0_0_0_yy_zz_0_yy = buffer_1000_ddsd[357];

    auto g_y_0_0_0_yy_zz_0_yz = buffer_1000_ddsd[358];

    auto g_y_0_0_0_yy_zz_0_zz = buffer_1000_ddsd[359];

    auto g_y_0_0_0_yz_xx_0_xx = buffer_1000_ddsd[360];

    auto g_y_0_0_0_yz_xx_0_xy = buffer_1000_ddsd[361];

    auto g_y_0_0_0_yz_xx_0_xz = buffer_1000_ddsd[362];

    auto g_y_0_0_0_yz_xx_0_yy = buffer_1000_ddsd[363];

    auto g_y_0_0_0_yz_xx_0_yz = buffer_1000_ddsd[364];

    auto g_y_0_0_0_yz_xx_0_zz = buffer_1000_ddsd[365];

    auto g_y_0_0_0_yz_xy_0_xx = buffer_1000_ddsd[366];

    auto g_y_0_0_0_yz_xy_0_xy = buffer_1000_ddsd[367];

    auto g_y_0_0_0_yz_xy_0_xz = buffer_1000_ddsd[368];

    auto g_y_0_0_0_yz_xy_0_yy = buffer_1000_ddsd[369];

    auto g_y_0_0_0_yz_xy_0_yz = buffer_1000_ddsd[370];

    auto g_y_0_0_0_yz_xy_0_zz = buffer_1000_ddsd[371];

    auto g_y_0_0_0_yz_xz_0_xx = buffer_1000_ddsd[372];

    auto g_y_0_0_0_yz_xz_0_xy = buffer_1000_ddsd[373];

    auto g_y_0_0_0_yz_xz_0_xz = buffer_1000_ddsd[374];

    auto g_y_0_0_0_yz_xz_0_yy = buffer_1000_ddsd[375];

    auto g_y_0_0_0_yz_xz_0_yz = buffer_1000_ddsd[376];

    auto g_y_0_0_0_yz_xz_0_zz = buffer_1000_ddsd[377];

    auto g_y_0_0_0_yz_yy_0_xx = buffer_1000_ddsd[378];

    auto g_y_0_0_0_yz_yy_0_xy = buffer_1000_ddsd[379];

    auto g_y_0_0_0_yz_yy_0_xz = buffer_1000_ddsd[380];

    auto g_y_0_0_0_yz_yy_0_yy = buffer_1000_ddsd[381];

    auto g_y_0_0_0_yz_yy_0_yz = buffer_1000_ddsd[382];

    auto g_y_0_0_0_yz_yy_0_zz = buffer_1000_ddsd[383];

    auto g_y_0_0_0_yz_yz_0_xx = buffer_1000_ddsd[384];

    auto g_y_0_0_0_yz_yz_0_xy = buffer_1000_ddsd[385];

    auto g_y_0_0_0_yz_yz_0_xz = buffer_1000_ddsd[386];

    auto g_y_0_0_0_yz_yz_0_yy = buffer_1000_ddsd[387];

    auto g_y_0_0_0_yz_yz_0_yz = buffer_1000_ddsd[388];

    auto g_y_0_0_0_yz_yz_0_zz = buffer_1000_ddsd[389];

    auto g_y_0_0_0_yz_zz_0_xx = buffer_1000_ddsd[390];

    auto g_y_0_0_0_yz_zz_0_xy = buffer_1000_ddsd[391];

    auto g_y_0_0_0_yz_zz_0_xz = buffer_1000_ddsd[392];

    auto g_y_0_0_0_yz_zz_0_yy = buffer_1000_ddsd[393];

    auto g_y_0_0_0_yz_zz_0_yz = buffer_1000_ddsd[394];

    auto g_y_0_0_0_yz_zz_0_zz = buffer_1000_ddsd[395];

    auto g_y_0_0_0_zz_xx_0_xx = buffer_1000_ddsd[396];

    auto g_y_0_0_0_zz_xx_0_xy = buffer_1000_ddsd[397];

    auto g_y_0_0_0_zz_xx_0_xz = buffer_1000_ddsd[398];

    auto g_y_0_0_0_zz_xx_0_yy = buffer_1000_ddsd[399];

    auto g_y_0_0_0_zz_xx_0_yz = buffer_1000_ddsd[400];

    auto g_y_0_0_0_zz_xx_0_zz = buffer_1000_ddsd[401];

    auto g_y_0_0_0_zz_xy_0_xx = buffer_1000_ddsd[402];

    auto g_y_0_0_0_zz_xy_0_xy = buffer_1000_ddsd[403];

    auto g_y_0_0_0_zz_xy_0_xz = buffer_1000_ddsd[404];

    auto g_y_0_0_0_zz_xy_0_yy = buffer_1000_ddsd[405];

    auto g_y_0_0_0_zz_xy_0_yz = buffer_1000_ddsd[406];

    auto g_y_0_0_0_zz_xy_0_zz = buffer_1000_ddsd[407];

    auto g_y_0_0_0_zz_xz_0_xx = buffer_1000_ddsd[408];

    auto g_y_0_0_0_zz_xz_0_xy = buffer_1000_ddsd[409];

    auto g_y_0_0_0_zz_xz_0_xz = buffer_1000_ddsd[410];

    auto g_y_0_0_0_zz_xz_0_yy = buffer_1000_ddsd[411];

    auto g_y_0_0_0_zz_xz_0_yz = buffer_1000_ddsd[412];

    auto g_y_0_0_0_zz_xz_0_zz = buffer_1000_ddsd[413];

    auto g_y_0_0_0_zz_yy_0_xx = buffer_1000_ddsd[414];

    auto g_y_0_0_0_zz_yy_0_xy = buffer_1000_ddsd[415];

    auto g_y_0_0_0_zz_yy_0_xz = buffer_1000_ddsd[416];

    auto g_y_0_0_0_zz_yy_0_yy = buffer_1000_ddsd[417];

    auto g_y_0_0_0_zz_yy_0_yz = buffer_1000_ddsd[418];

    auto g_y_0_0_0_zz_yy_0_zz = buffer_1000_ddsd[419];

    auto g_y_0_0_0_zz_yz_0_xx = buffer_1000_ddsd[420];

    auto g_y_0_0_0_zz_yz_0_xy = buffer_1000_ddsd[421];

    auto g_y_0_0_0_zz_yz_0_xz = buffer_1000_ddsd[422];

    auto g_y_0_0_0_zz_yz_0_yy = buffer_1000_ddsd[423];

    auto g_y_0_0_0_zz_yz_0_yz = buffer_1000_ddsd[424];

    auto g_y_0_0_0_zz_yz_0_zz = buffer_1000_ddsd[425];

    auto g_y_0_0_0_zz_zz_0_xx = buffer_1000_ddsd[426];

    auto g_y_0_0_0_zz_zz_0_xy = buffer_1000_ddsd[427];

    auto g_y_0_0_0_zz_zz_0_xz = buffer_1000_ddsd[428];

    auto g_y_0_0_0_zz_zz_0_yy = buffer_1000_ddsd[429];

    auto g_y_0_0_0_zz_zz_0_yz = buffer_1000_ddsd[430];

    auto g_y_0_0_0_zz_zz_0_zz = buffer_1000_ddsd[431];

    auto g_z_0_0_0_xx_xx_0_xx = buffer_1000_ddsd[432];

    auto g_z_0_0_0_xx_xx_0_xy = buffer_1000_ddsd[433];

    auto g_z_0_0_0_xx_xx_0_xz = buffer_1000_ddsd[434];

    auto g_z_0_0_0_xx_xx_0_yy = buffer_1000_ddsd[435];

    auto g_z_0_0_0_xx_xx_0_yz = buffer_1000_ddsd[436];

    auto g_z_0_0_0_xx_xx_0_zz = buffer_1000_ddsd[437];

    auto g_z_0_0_0_xx_xy_0_xx = buffer_1000_ddsd[438];

    auto g_z_0_0_0_xx_xy_0_xy = buffer_1000_ddsd[439];

    auto g_z_0_0_0_xx_xy_0_xz = buffer_1000_ddsd[440];

    auto g_z_0_0_0_xx_xy_0_yy = buffer_1000_ddsd[441];

    auto g_z_0_0_0_xx_xy_0_yz = buffer_1000_ddsd[442];

    auto g_z_0_0_0_xx_xy_0_zz = buffer_1000_ddsd[443];

    auto g_z_0_0_0_xx_xz_0_xx = buffer_1000_ddsd[444];

    auto g_z_0_0_0_xx_xz_0_xy = buffer_1000_ddsd[445];

    auto g_z_0_0_0_xx_xz_0_xz = buffer_1000_ddsd[446];

    auto g_z_0_0_0_xx_xz_0_yy = buffer_1000_ddsd[447];

    auto g_z_0_0_0_xx_xz_0_yz = buffer_1000_ddsd[448];

    auto g_z_0_0_0_xx_xz_0_zz = buffer_1000_ddsd[449];

    auto g_z_0_0_0_xx_yy_0_xx = buffer_1000_ddsd[450];

    auto g_z_0_0_0_xx_yy_0_xy = buffer_1000_ddsd[451];

    auto g_z_0_0_0_xx_yy_0_xz = buffer_1000_ddsd[452];

    auto g_z_0_0_0_xx_yy_0_yy = buffer_1000_ddsd[453];

    auto g_z_0_0_0_xx_yy_0_yz = buffer_1000_ddsd[454];

    auto g_z_0_0_0_xx_yy_0_zz = buffer_1000_ddsd[455];

    auto g_z_0_0_0_xx_yz_0_xx = buffer_1000_ddsd[456];

    auto g_z_0_0_0_xx_yz_0_xy = buffer_1000_ddsd[457];

    auto g_z_0_0_0_xx_yz_0_xz = buffer_1000_ddsd[458];

    auto g_z_0_0_0_xx_yz_0_yy = buffer_1000_ddsd[459];

    auto g_z_0_0_0_xx_yz_0_yz = buffer_1000_ddsd[460];

    auto g_z_0_0_0_xx_yz_0_zz = buffer_1000_ddsd[461];

    auto g_z_0_0_0_xx_zz_0_xx = buffer_1000_ddsd[462];

    auto g_z_0_0_0_xx_zz_0_xy = buffer_1000_ddsd[463];

    auto g_z_0_0_0_xx_zz_0_xz = buffer_1000_ddsd[464];

    auto g_z_0_0_0_xx_zz_0_yy = buffer_1000_ddsd[465];

    auto g_z_0_0_0_xx_zz_0_yz = buffer_1000_ddsd[466];

    auto g_z_0_0_0_xx_zz_0_zz = buffer_1000_ddsd[467];

    auto g_z_0_0_0_xy_xx_0_xx = buffer_1000_ddsd[468];

    auto g_z_0_0_0_xy_xx_0_xy = buffer_1000_ddsd[469];

    auto g_z_0_0_0_xy_xx_0_xz = buffer_1000_ddsd[470];

    auto g_z_0_0_0_xy_xx_0_yy = buffer_1000_ddsd[471];

    auto g_z_0_0_0_xy_xx_0_yz = buffer_1000_ddsd[472];

    auto g_z_0_0_0_xy_xx_0_zz = buffer_1000_ddsd[473];

    auto g_z_0_0_0_xy_xy_0_xx = buffer_1000_ddsd[474];

    auto g_z_0_0_0_xy_xy_0_xy = buffer_1000_ddsd[475];

    auto g_z_0_0_0_xy_xy_0_xz = buffer_1000_ddsd[476];

    auto g_z_0_0_0_xy_xy_0_yy = buffer_1000_ddsd[477];

    auto g_z_0_0_0_xy_xy_0_yz = buffer_1000_ddsd[478];

    auto g_z_0_0_0_xy_xy_0_zz = buffer_1000_ddsd[479];

    auto g_z_0_0_0_xy_xz_0_xx = buffer_1000_ddsd[480];

    auto g_z_0_0_0_xy_xz_0_xy = buffer_1000_ddsd[481];

    auto g_z_0_0_0_xy_xz_0_xz = buffer_1000_ddsd[482];

    auto g_z_0_0_0_xy_xz_0_yy = buffer_1000_ddsd[483];

    auto g_z_0_0_0_xy_xz_0_yz = buffer_1000_ddsd[484];

    auto g_z_0_0_0_xy_xz_0_zz = buffer_1000_ddsd[485];

    auto g_z_0_0_0_xy_yy_0_xx = buffer_1000_ddsd[486];

    auto g_z_0_0_0_xy_yy_0_xy = buffer_1000_ddsd[487];

    auto g_z_0_0_0_xy_yy_0_xz = buffer_1000_ddsd[488];

    auto g_z_0_0_0_xy_yy_0_yy = buffer_1000_ddsd[489];

    auto g_z_0_0_0_xy_yy_0_yz = buffer_1000_ddsd[490];

    auto g_z_0_0_0_xy_yy_0_zz = buffer_1000_ddsd[491];

    auto g_z_0_0_0_xy_yz_0_xx = buffer_1000_ddsd[492];

    auto g_z_0_0_0_xy_yz_0_xy = buffer_1000_ddsd[493];

    auto g_z_0_0_0_xy_yz_0_xz = buffer_1000_ddsd[494];

    auto g_z_0_0_0_xy_yz_0_yy = buffer_1000_ddsd[495];

    auto g_z_0_0_0_xy_yz_0_yz = buffer_1000_ddsd[496];

    auto g_z_0_0_0_xy_yz_0_zz = buffer_1000_ddsd[497];

    auto g_z_0_0_0_xy_zz_0_xx = buffer_1000_ddsd[498];

    auto g_z_0_0_0_xy_zz_0_xy = buffer_1000_ddsd[499];

    auto g_z_0_0_0_xy_zz_0_xz = buffer_1000_ddsd[500];

    auto g_z_0_0_0_xy_zz_0_yy = buffer_1000_ddsd[501];

    auto g_z_0_0_0_xy_zz_0_yz = buffer_1000_ddsd[502];

    auto g_z_0_0_0_xy_zz_0_zz = buffer_1000_ddsd[503];

    auto g_z_0_0_0_xz_xx_0_xx = buffer_1000_ddsd[504];

    auto g_z_0_0_0_xz_xx_0_xy = buffer_1000_ddsd[505];

    auto g_z_0_0_0_xz_xx_0_xz = buffer_1000_ddsd[506];

    auto g_z_0_0_0_xz_xx_0_yy = buffer_1000_ddsd[507];

    auto g_z_0_0_0_xz_xx_0_yz = buffer_1000_ddsd[508];

    auto g_z_0_0_0_xz_xx_0_zz = buffer_1000_ddsd[509];

    auto g_z_0_0_0_xz_xy_0_xx = buffer_1000_ddsd[510];

    auto g_z_0_0_0_xz_xy_0_xy = buffer_1000_ddsd[511];

    auto g_z_0_0_0_xz_xy_0_xz = buffer_1000_ddsd[512];

    auto g_z_0_0_0_xz_xy_0_yy = buffer_1000_ddsd[513];

    auto g_z_0_0_0_xz_xy_0_yz = buffer_1000_ddsd[514];

    auto g_z_0_0_0_xz_xy_0_zz = buffer_1000_ddsd[515];

    auto g_z_0_0_0_xz_xz_0_xx = buffer_1000_ddsd[516];

    auto g_z_0_0_0_xz_xz_0_xy = buffer_1000_ddsd[517];

    auto g_z_0_0_0_xz_xz_0_xz = buffer_1000_ddsd[518];

    auto g_z_0_0_0_xz_xz_0_yy = buffer_1000_ddsd[519];

    auto g_z_0_0_0_xz_xz_0_yz = buffer_1000_ddsd[520];

    auto g_z_0_0_0_xz_xz_0_zz = buffer_1000_ddsd[521];

    auto g_z_0_0_0_xz_yy_0_xx = buffer_1000_ddsd[522];

    auto g_z_0_0_0_xz_yy_0_xy = buffer_1000_ddsd[523];

    auto g_z_0_0_0_xz_yy_0_xz = buffer_1000_ddsd[524];

    auto g_z_0_0_0_xz_yy_0_yy = buffer_1000_ddsd[525];

    auto g_z_0_0_0_xz_yy_0_yz = buffer_1000_ddsd[526];

    auto g_z_0_0_0_xz_yy_0_zz = buffer_1000_ddsd[527];

    auto g_z_0_0_0_xz_yz_0_xx = buffer_1000_ddsd[528];

    auto g_z_0_0_0_xz_yz_0_xy = buffer_1000_ddsd[529];

    auto g_z_0_0_0_xz_yz_0_xz = buffer_1000_ddsd[530];

    auto g_z_0_0_0_xz_yz_0_yy = buffer_1000_ddsd[531];

    auto g_z_0_0_0_xz_yz_0_yz = buffer_1000_ddsd[532];

    auto g_z_0_0_0_xz_yz_0_zz = buffer_1000_ddsd[533];

    auto g_z_0_0_0_xz_zz_0_xx = buffer_1000_ddsd[534];

    auto g_z_0_0_0_xz_zz_0_xy = buffer_1000_ddsd[535];

    auto g_z_0_0_0_xz_zz_0_xz = buffer_1000_ddsd[536];

    auto g_z_0_0_0_xz_zz_0_yy = buffer_1000_ddsd[537];

    auto g_z_0_0_0_xz_zz_0_yz = buffer_1000_ddsd[538];

    auto g_z_0_0_0_xz_zz_0_zz = buffer_1000_ddsd[539];

    auto g_z_0_0_0_yy_xx_0_xx = buffer_1000_ddsd[540];

    auto g_z_0_0_0_yy_xx_0_xy = buffer_1000_ddsd[541];

    auto g_z_0_0_0_yy_xx_0_xz = buffer_1000_ddsd[542];

    auto g_z_0_0_0_yy_xx_0_yy = buffer_1000_ddsd[543];

    auto g_z_0_0_0_yy_xx_0_yz = buffer_1000_ddsd[544];

    auto g_z_0_0_0_yy_xx_0_zz = buffer_1000_ddsd[545];

    auto g_z_0_0_0_yy_xy_0_xx = buffer_1000_ddsd[546];

    auto g_z_0_0_0_yy_xy_0_xy = buffer_1000_ddsd[547];

    auto g_z_0_0_0_yy_xy_0_xz = buffer_1000_ddsd[548];

    auto g_z_0_0_0_yy_xy_0_yy = buffer_1000_ddsd[549];

    auto g_z_0_0_0_yy_xy_0_yz = buffer_1000_ddsd[550];

    auto g_z_0_0_0_yy_xy_0_zz = buffer_1000_ddsd[551];

    auto g_z_0_0_0_yy_xz_0_xx = buffer_1000_ddsd[552];

    auto g_z_0_0_0_yy_xz_0_xy = buffer_1000_ddsd[553];

    auto g_z_0_0_0_yy_xz_0_xz = buffer_1000_ddsd[554];

    auto g_z_0_0_0_yy_xz_0_yy = buffer_1000_ddsd[555];

    auto g_z_0_0_0_yy_xz_0_yz = buffer_1000_ddsd[556];

    auto g_z_0_0_0_yy_xz_0_zz = buffer_1000_ddsd[557];

    auto g_z_0_0_0_yy_yy_0_xx = buffer_1000_ddsd[558];

    auto g_z_0_0_0_yy_yy_0_xy = buffer_1000_ddsd[559];

    auto g_z_0_0_0_yy_yy_0_xz = buffer_1000_ddsd[560];

    auto g_z_0_0_0_yy_yy_0_yy = buffer_1000_ddsd[561];

    auto g_z_0_0_0_yy_yy_0_yz = buffer_1000_ddsd[562];

    auto g_z_0_0_0_yy_yy_0_zz = buffer_1000_ddsd[563];

    auto g_z_0_0_0_yy_yz_0_xx = buffer_1000_ddsd[564];

    auto g_z_0_0_0_yy_yz_0_xy = buffer_1000_ddsd[565];

    auto g_z_0_0_0_yy_yz_0_xz = buffer_1000_ddsd[566];

    auto g_z_0_0_0_yy_yz_0_yy = buffer_1000_ddsd[567];

    auto g_z_0_0_0_yy_yz_0_yz = buffer_1000_ddsd[568];

    auto g_z_0_0_0_yy_yz_0_zz = buffer_1000_ddsd[569];

    auto g_z_0_0_0_yy_zz_0_xx = buffer_1000_ddsd[570];

    auto g_z_0_0_0_yy_zz_0_xy = buffer_1000_ddsd[571];

    auto g_z_0_0_0_yy_zz_0_xz = buffer_1000_ddsd[572];

    auto g_z_0_0_0_yy_zz_0_yy = buffer_1000_ddsd[573];

    auto g_z_0_0_0_yy_zz_0_yz = buffer_1000_ddsd[574];

    auto g_z_0_0_0_yy_zz_0_zz = buffer_1000_ddsd[575];

    auto g_z_0_0_0_yz_xx_0_xx = buffer_1000_ddsd[576];

    auto g_z_0_0_0_yz_xx_0_xy = buffer_1000_ddsd[577];

    auto g_z_0_0_0_yz_xx_0_xz = buffer_1000_ddsd[578];

    auto g_z_0_0_0_yz_xx_0_yy = buffer_1000_ddsd[579];

    auto g_z_0_0_0_yz_xx_0_yz = buffer_1000_ddsd[580];

    auto g_z_0_0_0_yz_xx_0_zz = buffer_1000_ddsd[581];

    auto g_z_0_0_0_yz_xy_0_xx = buffer_1000_ddsd[582];

    auto g_z_0_0_0_yz_xy_0_xy = buffer_1000_ddsd[583];

    auto g_z_0_0_0_yz_xy_0_xz = buffer_1000_ddsd[584];

    auto g_z_0_0_0_yz_xy_0_yy = buffer_1000_ddsd[585];

    auto g_z_0_0_0_yz_xy_0_yz = buffer_1000_ddsd[586];

    auto g_z_0_0_0_yz_xy_0_zz = buffer_1000_ddsd[587];

    auto g_z_0_0_0_yz_xz_0_xx = buffer_1000_ddsd[588];

    auto g_z_0_0_0_yz_xz_0_xy = buffer_1000_ddsd[589];

    auto g_z_0_0_0_yz_xz_0_xz = buffer_1000_ddsd[590];

    auto g_z_0_0_0_yz_xz_0_yy = buffer_1000_ddsd[591];

    auto g_z_0_0_0_yz_xz_0_yz = buffer_1000_ddsd[592];

    auto g_z_0_0_0_yz_xz_0_zz = buffer_1000_ddsd[593];

    auto g_z_0_0_0_yz_yy_0_xx = buffer_1000_ddsd[594];

    auto g_z_0_0_0_yz_yy_0_xy = buffer_1000_ddsd[595];

    auto g_z_0_0_0_yz_yy_0_xz = buffer_1000_ddsd[596];

    auto g_z_0_0_0_yz_yy_0_yy = buffer_1000_ddsd[597];

    auto g_z_0_0_0_yz_yy_0_yz = buffer_1000_ddsd[598];

    auto g_z_0_0_0_yz_yy_0_zz = buffer_1000_ddsd[599];

    auto g_z_0_0_0_yz_yz_0_xx = buffer_1000_ddsd[600];

    auto g_z_0_0_0_yz_yz_0_xy = buffer_1000_ddsd[601];

    auto g_z_0_0_0_yz_yz_0_xz = buffer_1000_ddsd[602];

    auto g_z_0_0_0_yz_yz_0_yy = buffer_1000_ddsd[603];

    auto g_z_0_0_0_yz_yz_0_yz = buffer_1000_ddsd[604];

    auto g_z_0_0_0_yz_yz_0_zz = buffer_1000_ddsd[605];

    auto g_z_0_0_0_yz_zz_0_xx = buffer_1000_ddsd[606];

    auto g_z_0_0_0_yz_zz_0_xy = buffer_1000_ddsd[607];

    auto g_z_0_0_0_yz_zz_0_xz = buffer_1000_ddsd[608];

    auto g_z_0_0_0_yz_zz_0_yy = buffer_1000_ddsd[609];

    auto g_z_0_0_0_yz_zz_0_yz = buffer_1000_ddsd[610];

    auto g_z_0_0_0_yz_zz_0_zz = buffer_1000_ddsd[611];

    auto g_z_0_0_0_zz_xx_0_xx = buffer_1000_ddsd[612];

    auto g_z_0_0_0_zz_xx_0_xy = buffer_1000_ddsd[613];

    auto g_z_0_0_0_zz_xx_0_xz = buffer_1000_ddsd[614];

    auto g_z_0_0_0_zz_xx_0_yy = buffer_1000_ddsd[615];

    auto g_z_0_0_0_zz_xx_0_yz = buffer_1000_ddsd[616];

    auto g_z_0_0_0_zz_xx_0_zz = buffer_1000_ddsd[617];

    auto g_z_0_0_0_zz_xy_0_xx = buffer_1000_ddsd[618];

    auto g_z_0_0_0_zz_xy_0_xy = buffer_1000_ddsd[619];

    auto g_z_0_0_0_zz_xy_0_xz = buffer_1000_ddsd[620];

    auto g_z_0_0_0_zz_xy_0_yy = buffer_1000_ddsd[621];

    auto g_z_0_0_0_zz_xy_0_yz = buffer_1000_ddsd[622];

    auto g_z_0_0_0_zz_xy_0_zz = buffer_1000_ddsd[623];

    auto g_z_0_0_0_zz_xz_0_xx = buffer_1000_ddsd[624];

    auto g_z_0_0_0_zz_xz_0_xy = buffer_1000_ddsd[625];

    auto g_z_0_0_0_zz_xz_0_xz = buffer_1000_ddsd[626];

    auto g_z_0_0_0_zz_xz_0_yy = buffer_1000_ddsd[627];

    auto g_z_0_0_0_zz_xz_0_yz = buffer_1000_ddsd[628];

    auto g_z_0_0_0_zz_xz_0_zz = buffer_1000_ddsd[629];

    auto g_z_0_0_0_zz_yy_0_xx = buffer_1000_ddsd[630];

    auto g_z_0_0_0_zz_yy_0_xy = buffer_1000_ddsd[631];

    auto g_z_0_0_0_zz_yy_0_xz = buffer_1000_ddsd[632];

    auto g_z_0_0_0_zz_yy_0_yy = buffer_1000_ddsd[633];

    auto g_z_0_0_0_zz_yy_0_yz = buffer_1000_ddsd[634];

    auto g_z_0_0_0_zz_yy_0_zz = buffer_1000_ddsd[635];

    auto g_z_0_0_0_zz_yz_0_xx = buffer_1000_ddsd[636];

    auto g_z_0_0_0_zz_yz_0_xy = buffer_1000_ddsd[637];

    auto g_z_0_0_0_zz_yz_0_xz = buffer_1000_ddsd[638];

    auto g_z_0_0_0_zz_yz_0_yy = buffer_1000_ddsd[639];

    auto g_z_0_0_0_zz_yz_0_yz = buffer_1000_ddsd[640];

    auto g_z_0_0_0_zz_yz_0_zz = buffer_1000_ddsd[641];

    auto g_z_0_0_0_zz_zz_0_xx = buffer_1000_ddsd[642];

    auto g_z_0_0_0_zz_zz_0_xy = buffer_1000_ddsd[643];

    auto g_z_0_0_0_zz_zz_0_xz = buffer_1000_ddsd[644];

    auto g_z_0_0_0_zz_zz_0_yy = buffer_1000_ddsd[645];

    auto g_z_0_0_0_zz_zz_0_yz = buffer_1000_ddsd[646];

    auto g_z_0_0_0_zz_zz_0_zz = buffer_1000_ddsd[647];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_0_xx, g_x_0_0_0_xx_xx_0_xy, g_x_0_0_0_xx_xx_0_xz, g_x_0_0_0_xx_xx_0_yy, g_x_0_0_0_xx_xx_0_yz, g_x_0_0_0_xx_xx_0_zz, g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_xxx_xx_0_xx, g_xxx_xx_0_xy, g_xxx_xx_0_xz, g_xxx_xx_0_yy, g_xxx_xx_0_yz, g_xxx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_0_xx[i] = -2.0 * g_x_xx_0_xx[i] + 2.0 * g_xxx_xx_0_xx[i] * a_exp;

        g_x_0_0_0_xx_xx_0_xy[i] = -2.0 * g_x_xx_0_xy[i] + 2.0 * g_xxx_xx_0_xy[i] * a_exp;

        g_x_0_0_0_xx_xx_0_xz[i] = -2.0 * g_x_xx_0_xz[i] + 2.0 * g_xxx_xx_0_xz[i] * a_exp;

        g_x_0_0_0_xx_xx_0_yy[i] = -2.0 * g_x_xx_0_yy[i] + 2.0 * g_xxx_xx_0_yy[i] * a_exp;

        g_x_0_0_0_xx_xx_0_yz[i] = -2.0 * g_x_xx_0_yz[i] + 2.0 * g_xxx_xx_0_yz[i] * a_exp;

        g_x_0_0_0_xx_xx_0_zz[i] = -2.0 * g_x_xx_0_zz[i] + 2.0 * g_xxx_xx_0_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_0_xx, g_x_0_0_0_xx_xy_0_xy, g_x_0_0_0_xx_xy_0_xz, g_x_0_0_0_xx_xy_0_yy, g_x_0_0_0_xx_xy_0_yz, g_x_0_0_0_xx_xy_0_zz, g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_xxx_xy_0_xx, g_xxx_xy_0_xy, g_xxx_xy_0_xz, g_xxx_xy_0_yy, g_xxx_xy_0_yz, g_xxx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_0_xx[i] = -2.0 * g_x_xy_0_xx[i] + 2.0 * g_xxx_xy_0_xx[i] * a_exp;

        g_x_0_0_0_xx_xy_0_xy[i] = -2.0 * g_x_xy_0_xy[i] + 2.0 * g_xxx_xy_0_xy[i] * a_exp;

        g_x_0_0_0_xx_xy_0_xz[i] = -2.0 * g_x_xy_0_xz[i] + 2.0 * g_xxx_xy_0_xz[i] * a_exp;

        g_x_0_0_0_xx_xy_0_yy[i] = -2.0 * g_x_xy_0_yy[i] + 2.0 * g_xxx_xy_0_yy[i] * a_exp;

        g_x_0_0_0_xx_xy_0_yz[i] = -2.0 * g_x_xy_0_yz[i] + 2.0 * g_xxx_xy_0_yz[i] * a_exp;

        g_x_0_0_0_xx_xy_0_zz[i] = -2.0 * g_x_xy_0_zz[i] + 2.0 * g_xxx_xy_0_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_0_xx, g_x_0_0_0_xx_xz_0_xy, g_x_0_0_0_xx_xz_0_xz, g_x_0_0_0_xx_xz_0_yy, g_x_0_0_0_xx_xz_0_yz, g_x_0_0_0_xx_xz_0_zz, g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_xxx_xz_0_xx, g_xxx_xz_0_xy, g_xxx_xz_0_xz, g_xxx_xz_0_yy, g_xxx_xz_0_yz, g_xxx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_0_xx[i] = -2.0 * g_x_xz_0_xx[i] + 2.0 * g_xxx_xz_0_xx[i] * a_exp;

        g_x_0_0_0_xx_xz_0_xy[i] = -2.0 * g_x_xz_0_xy[i] + 2.0 * g_xxx_xz_0_xy[i] * a_exp;

        g_x_0_0_0_xx_xz_0_xz[i] = -2.0 * g_x_xz_0_xz[i] + 2.0 * g_xxx_xz_0_xz[i] * a_exp;

        g_x_0_0_0_xx_xz_0_yy[i] = -2.0 * g_x_xz_0_yy[i] + 2.0 * g_xxx_xz_0_yy[i] * a_exp;

        g_x_0_0_0_xx_xz_0_yz[i] = -2.0 * g_x_xz_0_yz[i] + 2.0 * g_xxx_xz_0_yz[i] * a_exp;

        g_x_0_0_0_xx_xz_0_zz[i] = -2.0 * g_x_xz_0_zz[i] + 2.0 * g_xxx_xz_0_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_0_xx, g_x_0_0_0_xx_yy_0_xy, g_x_0_0_0_xx_yy_0_xz, g_x_0_0_0_xx_yy_0_yy, g_x_0_0_0_xx_yy_0_yz, g_x_0_0_0_xx_yy_0_zz, g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_xxx_yy_0_xx, g_xxx_yy_0_xy, g_xxx_yy_0_xz, g_xxx_yy_0_yy, g_xxx_yy_0_yz, g_xxx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_0_xx[i] = -2.0 * g_x_yy_0_xx[i] + 2.0 * g_xxx_yy_0_xx[i] * a_exp;

        g_x_0_0_0_xx_yy_0_xy[i] = -2.0 * g_x_yy_0_xy[i] + 2.0 * g_xxx_yy_0_xy[i] * a_exp;

        g_x_0_0_0_xx_yy_0_xz[i] = -2.0 * g_x_yy_0_xz[i] + 2.0 * g_xxx_yy_0_xz[i] * a_exp;

        g_x_0_0_0_xx_yy_0_yy[i] = -2.0 * g_x_yy_0_yy[i] + 2.0 * g_xxx_yy_0_yy[i] * a_exp;

        g_x_0_0_0_xx_yy_0_yz[i] = -2.0 * g_x_yy_0_yz[i] + 2.0 * g_xxx_yy_0_yz[i] * a_exp;

        g_x_0_0_0_xx_yy_0_zz[i] = -2.0 * g_x_yy_0_zz[i] + 2.0 * g_xxx_yy_0_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_0_xx, g_x_0_0_0_xx_yz_0_xy, g_x_0_0_0_xx_yz_0_xz, g_x_0_0_0_xx_yz_0_yy, g_x_0_0_0_xx_yz_0_yz, g_x_0_0_0_xx_yz_0_zz, g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_xxx_yz_0_xx, g_xxx_yz_0_xy, g_xxx_yz_0_xz, g_xxx_yz_0_yy, g_xxx_yz_0_yz, g_xxx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_0_xx[i] = -2.0 * g_x_yz_0_xx[i] + 2.0 * g_xxx_yz_0_xx[i] * a_exp;

        g_x_0_0_0_xx_yz_0_xy[i] = -2.0 * g_x_yz_0_xy[i] + 2.0 * g_xxx_yz_0_xy[i] * a_exp;

        g_x_0_0_0_xx_yz_0_xz[i] = -2.0 * g_x_yz_0_xz[i] + 2.0 * g_xxx_yz_0_xz[i] * a_exp;

        g_x_0_0_0_xx_yz_0_yy[i] = -2.0 * g_x_yz_0_yy[i] + 2.0 * g_xxx_yz_0_yy[i] * a_exp;

        g_x_0_0_0_xx_yz_0_yz[i] = -2.0 * g_x_yz_0_yz[i] + 2.0 * g_xxx_yz_0_yz[i] * a_exp;

        g_x_0_0_0_xx_yz_0_zz[i] = -2.0 * g_x_yz_0_zz[i] + 2.0 * g_xxx_yz_0_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_0_xx, g_x_0_0_0_xx_zz_0_xy, g_x_0_0_0_xx_zz_0_xz, g_x_0_0_0_xx_zz_0_yy, g_x_0_0_0_xx_zz_0_yz, g_x_0_0_0_xx_zz_0_zz, g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_xxx_zz_0_xx, g_xxx_zz_0_xy, g_xxx_zz_0_xz, g_xxx_zz_0_yy, g_xxx_zz_0_yz, g_xxx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_0_xx[i] = -2.0 * g_x_zz_0_xx[i] + 2.0 * g_xxx_zz_0_xx[i] * a_exp;

        g_x_0_0_0_xx_zz_0_xy[i] = -2.0 * g_x_zz_0_xy[i] + 2.0 * g_xxx_zz_0_xy[i] * a_exp;

        g_x_0_0_0_xx_zz_0_xz[i] = -2.0 * g_x_zz_0_xz[i] + 2.0 * g_xxx_zz_0_xz[i] * a_exp;

        g_x_0_0_0_xx_zz_0_yy[i] = -2.0 * g_x_zz_0_yy[i] + 2.0 * g_xxx_zz_0_yy[i] * a_exp;

        g_x_0_0_0_xx_zz_0_yz[i] = -2.0 * g_x_zz_0_yz[i] + 2.0 * g_xxx_zz_0_yz[i] * a_exp;

        g_x_0_0_0_xx_zz_0_zz[i] = -2.0 * g_x_zz_0_zz[i] + 2.0 * g_xxx_zz_0_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_0_xx, g_x_0_0_0_xy_xx_0_xy, g_x_0_0_0_xy_xx_0_xz, g_x_0_0_0_xy_xx_0_yy, g_x_0_0_0_xy_xx_0_yz, g_x_0_0_0_xy_xx_0_zz, g_xxy_xx_0_xx, g_xxy_xx_0_xy, g_xxy_xx_0_xz, g_xxy_xx_0_yy, g_xxy_xx_0_yz, g_xxy_xx_0_zz, g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_0_xx[i] = -g_y_xx_0_xx[i] + 2.0 * g_xxy_xx_0_xx[i] * a_exp;

        g_x_0_0_0_xy_xx_0_xy[i] = -g_y_xx_0_xy[i] + 2.0 * g_xxy_xx_0_xy[i] * a_exp;

        g_x_0_0_0_xy_xx_0_xz[i] = -g_y_xx_0_xz[i] + 2.0 * g_xxy_xx_0_xz[i] * a_exp;

        g_x_0_0_0_xy_xx_0_yy[i] = -g_y_xx_0_yy[i] + 2.0 * g_xxy_xx_0_yy[i] * a_exp;

        g_x_0_0_0_xy_xx_0_yz[i] = -g_y_xx_0_yz[i] + 2.0 * g_xxy_xx_0_yz[i] * a_exp;

        g_x_0_0_0_xy_xx_0_zz[i] = -g_y_xx_0_zz[i] + 2.0 * g_xxy_xx_0_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_0_xx, g_x_0_0_0_xy_xy_0_xy, g_x_0_0_0_xy_xy_0_xz, g_x_0_0_0_xy_xy_0_yy, g_x_0_0_0_xy_xy_0_yz, g_x_0_0_0_xy_xy_0_zz, g_xxy_xy_0_xx, g_xxy_xy_0_xy, g_xxy_xy_0_xz, g_xxy_xy_0_yy, g_xxy_xy_0_yz, g_xxy_xy_0_zz, g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_0_xx[i] = -g_y_xy_0_xx[i] + 2.0 * g_xxy_xy_0_xx[i] * a_exp;

        g_x_0_0_0_xy_xy_0_xy[i] = -g_y_xy_0_xy[i] + 2.0 * g_xxy_xy_0_xy[i] * a_exp;

        g_x_0_0_0_xy_xy_0_xz[i] = -g_y_xy_0_xz[i] + 2.0 * g_xxy_xy_0_xz[i] * a_exp;

        g_x_0_0_0_xy_xy_0_yy[i] = -g_y_xy_0_yy[i] + 2.0 * g_xxy_xy_0_yy[i] * a_exp;

        g_x_0_0_0_xy_xy_0_yz[i] = -g_y_xy_0_yz[i] + 2.0 * g_xxy_xy_0_yz[i] * a_exp;

        g_x_0_0_0_xy_xy_0_zz[i] = -g_y_xy_0_zz[i] + 2.0 * g_xxy_xy_0_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_0_xx, g_x_0_0_0_xy_xz_0_xy, g_x_0_0_0_xy_xz_0_xz, g_x_0_0_0_xy_xz_0_yy, g_x_0_0_0_xy_xz_0_yz, g_x_0_0_0_xy_xz_0_zz, g_xxy_xz_0_xx, g_xxy_xz_0_xy, g_xxy_xz_0_xz, g_xxy_xz_0_yy, g_xxy_xz_0_yz, g_xxy_xz_0_zz, g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_0_xx[i] = -g_y_xz_0_xx[i] + 2.0 * g_xxy_xz_0_xx[i] * a_exp;

        g_x_0_0_0_xy_xz_0_xy[i] = -g_y_xz_0_xy[i] + 2.0 * g_xxy_xz_0_xy[i] * a_exp;

        g_x_0_0_0_xy_xz_0_xz[i] = -g_y_xz_0_xz[i] + 2.0 * g_xxy_xz_0_xz[i] * a_exp;

        g_x_0_0_0_xy_xz_0_yy[i] = -g_y_xz_0_yy[i] + 2.0 * g_xxy_xz_0_yy[i] * a_exp;

        g_x_0_0_0_xy_xz_0_yz[i] = -g_y_xz_0_yz[i] + 2.0 * g_xxy_xz_0_yz[i] * a_exp;

        g_x_0_0_0_xy_xz_0_zz[i] = -g_y_xz_0_zz[i] + 2.0 * g_xxy_xz_0_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_0_xx, g_x_0_0_0_xy_yy_0_xy, g_x_0_0_0_xy_yy_0_xz, g_x_0_0_0_xy_yy_0_yy, g_x_0_0_0_xy_yy_0_yz, g_x_0_0_0_xy_yy_0_zz, g_xxy_yy_0_xx, g_xxy_yy_0_xy, g_xxy_yy_0_xz, g_xxy_yy_0_yy, g_xxy_yy_0_yz, g_xxy_yy_0_zz, g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_0_xx[i] = -g_y_yy_0_xx[i] + 2.0 * g_xxy_yy_0_xx[i] * a_exp;

        g_x_0_0_0_xy_yy_0_xy[i] = -g_y_yy_0_xy[i] + 2.0 * g_xxy_yy_0_xy[i] * a_exp;

        g_x_0_0_0_xy_yy_0_xz[i] = -g_y_yy_0_xz[i] + 2.0 * g_xxy_yy_0_xz[i] * a_exp;

        g_x_0_0_0_xy_yy_0_yy[i] = -g_y_yy_0_yy[i] + 2.0 * g_xxy_yy_0_yy[i] * a_exp;

        g_x_0_0_0_xy_yy_0_yz[i] = -g_y_yy_0_yz[i] + 2.0 * g_xxy_yy_0_yz[i] * a_exp;

        g_x_0_0_0_xy_yy_0_zz[i] = -g_y_yy_0_zz[i] + 2.0 * g_xxy_yy_0_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_0_xx, g_x_0_0_0_xy_yz_0_xy, g_x_0_0_0_xy_yz_0_xz, g_x_0_0_0_xy_yz_0_yy, g_x_0_0_0_xy_yz_0_yz, g_x_0_0_0_xy_yz_0_zz, g_xxy_yz_0_xx, g_xxy_yz_0_xy, g_xxy_yz_0_xz, g_xxy_yz_0_yy, g_xxy_yz_0_yz, g_xxy_yz_0_zz, g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_0_xx[i] = -g_y_yz_0_xx[i] + 2.0 * g_xxy_yz_0_xx[i] * a_exp;

        g_x_0_0_0_xy_yz_0_xy[i] = -g_y_yz_0_xy[i] + 2.0 * g_xxy_yz_0_xy[i] * a_exp;

        g_x_0_0_0_xy_yz_0_xz[i] = -g_y_yz_0_xz[i] + 2.0 * g_xxy_yz_0_xz[i] * a_exp;

        g_x_0_0_0_xy_yz_0_yy[i] = -g_y_yz_0_yy[i] + 2.0 * g_xxy_yz_0_yy[i] * a_exp;

        g_x_0_0_0_xy_yz_0_yz[i] = -g_y_yz_0_yz[i] + 2.0 * g_xxy_yz_0_yz[i] * a_exp;

        g_x_0_0_0_xy_yz_0_zz[i] = -g_y_yz_0_zz[i] + 2.0 * g_xxy_yz_0_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_0_xx, g_x_0_0_0_xy_zz_0_xy, g_x_0_0_0_xy_zz_0_xz, g_x_0_0_0_xy_zz_0_yy, g_x_0_0_0_xy_zz_0_yz, g_x_0_0_0_xy_zz_0_zz, g_xxy_zz_0_xx, g_xxy_zz_0_xy, g_xxy_zz_0_xz, g_xxy_zz_0_yy, g_xxy_zz_0_yz, g_xxy_zz_0_zz, g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_0_xx[i] = -g_y_zz_0_xx[i] + 2.0 * g_xxy_zz_0_xx[i] * a_exp;

        g_x_0_0_0_xy_zz_0_xy[i] = -g_y_zz_0_xy[i] + 2.0 * g_xxy_zz_0_xy[i] * a_exp;

        g_x_0_0_0_xy_zz_0_xz[i] = -g_y_zz_0_xz[i] + 2.0 * g_xxy_zz_0_xz[i] * a_exp;

        g_x_0_0_0_xy_zz_0_yy[i] = -g_y_zz_0_yy[i] + 2.0 * g_xxy_zz_0_yy[i] * a_exp;

        g_x_0_0_0_xy_zz_0_yz[i] = -g_y_zz_0_yz[i] + 2.0 * g_xxy_zz_0_yz[i] * a_exp;

        g_x_0_0_0_xy_zz_0_zz[i] = -g_y_zz_0_zz[i] + 2.0 * g_xxy_zz_0_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_0_xx, g_x_0_0_0_xz_xx_0_xy, g_x_0_0_0_xz_xx_0_xz, g_x_0_0_0_xz_xx_0_yy, g_x_0_0_0_xz_xx_0_yz, g_x_0_0_0_xz_xx_0_zz, g_xxz_xx_0_xx, g_xxz_xx_0_xy, g_xxz_xx_0_xz, g_xxz_xx_0_yy, g_xxz_xx_0_yz, g_xxz_xx_0_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_0_xx[i] = -g_z_xx_0_xx[i] + 2.0 * g_xxz_xx_0_xx[i] * a_exp;

        g_x_0_0_0_xz_xx_0_xy[i] = -g_z_xx_0_xy[i] + 2.0 * g_xxz_xx_0_xy[i] * a_exp;

        g_x_0_0_0_xz_xx_0_xz[i] = -g_z_xx_0_xz[i] + 2.0 * g_xxz_xx_0_xz[i] * a_exp;

        g_x_0_0_0_xz_xx_0_yy[i] = -g_z_xx_0_yy[i] + 2.0 * g_xxz_xx_0_yy[i] * a_exp;

        g_x_0_0_0_xz_xx_0_yz[i] = -g_z_xx_0_yz[i] + 2.0 * g_xxz_xx_0_yz[i] * a_exp;

        g_x_0_0_0_xz_xx_0_zz[i] = -g_z_xx_0_zz[i] + 2.0 * g_xxz_xx_0_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_0_xx, g_x_0_0_0_xz_xy_0_xy, g_x_0_0_0_xz_xy_0_xz, g_x_0_0_0_xz_xy_0_yy, g_x_0_0_0_xz_xy_0_yz, g_x_0_0_0_xz_xy_0_zz, g_xxz_xy_0_xx, g_xxz_xy_0_xy, g_xxz_xy_0_xz, g_xxz_xy_0_yy, g_xxz_xy_0_yz, g_xxz_xy_0_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_0_xx[i] = -g_z_xy_0_xx[i] + 2.0 * g_xxz_xy_0_xx[i] * a_exp;

        g_x_0_0_0_xz_xy_0_xy[i] = -g_z_xy_0_xy[i] + 2.0 * g_xxz_xy_0_xy[i] * a_exp;

        g_x_0_0_0_xz_xy_0_xz[i] = -g_z_xy_0_xz[i] + 2.0 * g_xxz_xy_0_xz[i] * a_exp;

        g_x_0_0_0_xz_xy_0_yy[i] = -g_z_xy_0_yy[i] + 2.0 * g_xxz_xy_0_yy[i] * a_exp;

        g_x_0_0_0_xz_xy_0_yz[i] = -g_z_xy_0_yz[i] + 2.0 * g_xxz_xy_0_yz[i] * a_exp;

        g_x_0_0_0_xz_xy_0_zz[i] = -g_z_xy_0_zz[i] + 2.0 * g_xxz_xy_0_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_0_xx, g_x_0_0_0_xz_xz_0_xy, g_x_0_0_0_xz_xz_0_xz, g_x_0_0_0_xz_xz_0_yy, g_x_0_0_0_xz_xz_0_yz, g_x_0_0_0_xz_xz_0_zz, g_xxz_xz_0_xx, g_xxz_xz_0_xy, g_xxz_xz_0_xz, g_xxz_xz_0_yy, g_xxz_xz_0_yz, g_xxz_xz_0_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_0_xx[i] = -g_z_xz_0_xx[i] + 2.0 * g_xxz_xz_0_xx[i] * a_exp;

        g_x_0_0_0_xz_xz_0_xy[i] = -g_z_xz_0_xy[i] + 2.0 * g_xxz_xz_0_xy[i] * a_exp;

        g_x_0_0_0_xz_xz_0_xz[i] = -g_z_xz_0_xz[i] + 2.0 * g_xxz_xz_0_xz[i] * a_exp;

        g_x_0_0_0_xz_xz_0_yy[i] = -g_z_xz_0_yy[i] + 2.0 * g_xxz_xz_0_yy[i] * a_exp;

        g_x_0_0_0_xz_xz_0_yz[i] = -g_z_xz_0_yz[i] + 2.0 * g_xxz_xz_0_yz[i] * a_exp;

        g_x_0_0_0_xz_xz_0_zz[i] = -g_z_xz_0_zz[i] + 2.0 * g_xxz_xz_0_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_0_xx, g_x_0_0_0_xz_yy_0_xy, g_x_0_0_0_xz_yy_0_xz, g_x_0_0_0_xz_yy_0_yy, g_x_0_0_0_xz_yy_0_yz, g_x_0_0_0_xz_yy_0_zz, g_xxz_yy_0_xx, g_xxz_yy_0_xy, g_xxz_yy_0_xz, g_xxz_yy_0_yy, g_xxz_yy_0_yz, g_xxz_yy_0_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_0_xx[i] = -g_z_yy_0_xx[i] + 2.0 * g_xxz_yy_0_xx[i] * a_exp;

        g_x_0_0_0_xz_yy_0_xy[i] = -g_z_yy_0_xy[i] + 2.0 * g_xxz_yy_0_xy[i] * a_exp;

        g_x_0_0_0_xz_yy_0_xz[i] = -g_z_yy_0_xz[i] + 2.0 * g_xxz_yy_0_xz[i] * a_exp;

        g_x_0_0_0_xz_yy_0_yy[i] = -g_z_yy_0_yy[i] + 2.0 * g_xxz_yy_0_yy[i] * a_exp;

        g_x_0_0_0_xz_yy_0_yz[i] = -g_z_yy_0_yz[i] + 2.0 * g_xxz_yy_0_yz[i] * a_exp;

        g_x_0_0_0_xz_yy_0_zz[i] = -g_z_yy_0_zz[i] + 2.0 * g_xxz_yy_0_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_0_xx, g_x_0_0_0_xz_yz_0_xy, g_x_0_0_0_xz_yz_0_xz, g_x_0_0_0_xz_yz_0_yy, g_x_0_0_0_xz_yz_0_yz, g_x_0_0_0_xz_yz_0_zz, g_xxz_yz_0_xx, g_xxz_yz_0_xy, g_xxz_yz_0_xz, g_xxz_yz_0_yy, g_xxz_yz_0_yz, g_xxz_yz_0_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_0_xx[i] = -g_z_yz_0_xx[i] + 2.0 * g_xxz_yz_0_xx[i] * a_exp;

        g_x_0_0_0_xz_yz_0_xy[i] = -g_z_yz_0_xy[i] + 2.0 * g_xxz_yz_0_xy[i] * a_exp;

        g_x_0_0_0_xz_yz_0_xz[i] = -g_z_yz_0_xz[i] + 2.0 * g_xxz_yz_0_xz[i] * a_exp;

        g_x_0_0_0_xz_yz_0_yy[i] = -g_z_yz_0_yy[i] + 2.0 * g_xxz_yz_0_yy[i] * a_exp;

        g_x_0_0_0_xz_yz_0_yz[i] = -g_z_yz_0_yz[i] + 2.0 * g_xxz_yz_0_yz[i] * a_exp;

        g_x_0_0_0_xz_yz_0_zz[i] = -g_z_yz_0_zz[i] + 2.0 * g_xxz_yz_0_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_0_xx, g_x_0_0_0_xz_zz_0_xy, g_x_0_0_0_xz_zz_0_xz, g_x_0_0_0_xz_zz_0_yy, g_x_0_0_0_xz_zz_0_yz, g_x_0_0_0_xz_zz_0_zz, g_xxz_zz_0_xx, g_xxz_zz_0_xy, g_xxz_zz_0_xz, g_xxz_zz_0_yy, g_xxz_zz_0_yz, g_xxz_zz_0_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_0_xx[i] = -g_z_zz_0_xx[i] + 2.0 * g_xxz_zz_0_xx[i] * a_exp;

        g_x_0_0_0_xz_zz_0_xy[i] = -g_z_zz_0_xy[i] + 2.0 * g_xxz_zz_0_xy[i] * a_exp;

        g_x_0_0_0_xz_zz_0_xz[i] = -g_z_zz_0_xz[i] + 2.0 * g_xxz_zz_0_xz[i] * a_exp;

        g_x_0_0_0_xz_zz_0_yy[i] = -g_z_zz_0_yy[i] + 2.0 * g_xxz_zz_0_yy[i] * a_exp;

        g_x_0_0_0_xz_zz_0_yz[i] = -g_z_zz_0_yz[i] + 2.0 * g_xxz_zz_0_yz[i] * a_exp;

        g_x_0_0_0_xz_zz_0_zz[i] = -g_z_zz_0_zz[i] + 2.0 * g_xxz_zz_0_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_0_xx, g_x_0_0_0_yy_xx_0_xy, g_x_0_0_0_yy_xx_0_xz, g_x_0_0_0_yy_xx_0_yy, g_x_0_0_0_yy_xx_0_yz, g_x_0_0_0_yy_xx_0_zz, g_xyy_xx_0_xx, g_xyy_xx_0_xy, g_xyy_xx_0_xz, g_xyy_xx_0_yy, g_xyy_xx_0_yz, g_xyy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_0_xx[i] = 2.0 * g_xyy_xx_0_xx[i] * a_exp;

        g_x_0_0_0_yy_xx_0_xy[i] = 2.0 * g_xyy_xx_0_xy[i] * a_exp;

        g_x_0_0_0_yy_xx_0_xz[i] = 2.0 * g_xyy_xx_0_xz[i] * a_exp;

        g_x_0_0_0_yy_xx_0_yy[i] = 2.0 * g_xyy_xx_0_yy[i] * a_exp;

        g_x_0_0_0_yy_xx_0_yz[i] = 2.0 * g_xyy_xx_0_yz[i] * a_exp;

        g_x_0_0_0_yy_xx_0_zz[i] = 2.0 * g_xyy_xx_0_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_0_xx, g_x_0_0_0_yy_xy_0_xy, g_x_0_0_0_yy_xy_0_xz, g_x_0_0_0_yy_xy_0_yy, g_x_0_0_0_yy_xy_0_yz, g_x_0_0_0_yy_xy_0_zz, g_xyy_xy_0_xx, g_xyy_xy_0_xy, g_xyy_xy_0_xz, g_xyy_xy_0_yy, g_xyy_xy_0_yz, g_xyy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_0_xx[i] = 2.0 * g_xyy_xy_0_xx[i] * a_exp;

        g_x_0_0_0_yy_xy_0_xy[i] = 2.0 * g_xyy_xy_0_xy[i] * a_exp;

        g_x_0_0_0_yy_xy_0_xz[i] = 2.0 * g_xyy_xy_0_xz[i] * a_exp;

        g_x_0_0_0_yy_xy_0_yy[i] = 2.0 * g_xyy_xy_0_yy[i] * a_exp;

        g_x_0_0_0_yy_xy_0_yz[i] = 2.0 * g_xyy_xy_0_yz[i] * a_exp;

        g_x_0_0_0_yy_xy_0_zz[i] = 2.0 * g_xyy_xy_0_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_0_xx, g_x_0_0_0_yy_xz_0_xy, g_x_0_0_0_yy_xz_0_xz, g_x_0_0_0_yy_xz_0_yy, g_x_0_0_0_yy_xz_0_yz, g_x_0_0_0_yy_xz_0_zz, g_xyy_xz_0_xx, g_xyy_xz_0_xy, g_xyy_xz_0_xz, g_xyy_xz_0_yy, g_xyy_xz_0_yz, g_xyy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_0_xx[i] = 2.0 * g_xyy_xz_0_xx[i] * a_exp;

        g_x_0_0_0_yy_xz_0_xy[i] = 2.0 * g_xyy_xz_0_xy[i] * a_exp;

        g_x_0_0_0_yy_xz_0_xz[i] = 2.0 * g_xyy_xz_0_xz[i] * a_exp;

        g_x_0_0_0_yy_xz_0_yy[i] = 2.0 * g_xyy_xz_0_yy[i] * a_exp;

        g_x_0_0_0_yy_xz_0_yz[i] = 2.0 * g_xyy_xz_0_yz[i] * a_exp;

        g_x_0_0_0_yy_xz_0_zz[i] = 2.0 * g_xyy_xz_0_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_0_xx, g_x_0_0_0_yy_yy_0_xy, g_x_0_0_0_yy_yy_0_xz, g_x_0_0_0_yy_yy_0_yy, g_x_0_0_0_yy_yy_0_yz, g_x_0_0_0_yy_yy_0_zz, g_xyy_yy_0_xx, g_xyy_yy_0_xy, g_xyy_yy_0_xz, g_xyy_yy_0_yy, g_xyy_yy_0_yz, g_xyy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_0_xx[i] = 2.0 * g_xyy_yy_0_xx[i] * a_exp;

        g_x_0_0_0_yy_yy_0_xy[i] = 2.0 * g_xyy_yy_0_xy[i] * a_exp;

        g_x_0_0_0_yy_yy_0_xz[i] = 2.0 * g_xyy_yy_0_xz[i] * a_exp;

        g_x_0_0_0_yy_yy_0_yy[i] = 2.0 * g_xyy_yy_0_yy[i] * a_exp;

        g_x_0_0_0_yy_yy_0_yz[i] = 2.0 * g_xyy_yy_0_yz[i] * a_exp;

        g_x_0_0_0_yy_yy_0_zz[i] = 2.0 * g_xyy_yy_0_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_0_xx, g_x_0_0_0_yy_yz_0_xy, g_x_0_0_0_yy_yz_0_xz, g_x_0_0_0_yy_yz_0_yy, g_x_0_0_0_yy_yz_0_yz, g_x_0_0_0_yy_yz_0_zz, g_xyy_yz_0_xx, g_xyy_yz_0_xy, g_xyy_yz_0_xz, g_xyy_yz_0_yy, g_xyy_yz_0_yz, g_xyy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_0_xx[i] = 2.0 * g_xyy_yz_0_xx[i] * a_exp;

        g_x_0_0_0_yy_yz_0_xy[i] = 2.0 * g_xyy_yz_0_xy[i] * a_exp;

        g_x_0_0_0_yy_yz_0_xz[i] = 2.0 * g_xyy_yz_0_xz[i] * a_exp;

        g_x_0_0_0_yy_yz_0_yy[i] = 2.0 * g_xyy_yz_0_yy[i] * a_exp;

        g_x_0_0_0_yy_yz_0_yz[i] = 2.0 * g_xyy_yz_0_yz[i] * a_exp;

        g_x_0_0_0_yy_yz_0_zz[i] = 2.0 * g_xyy_yz_0_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_0_xx, g_x_0_0_0_yy_zz_0_xy, g_x_0_0_0_yy_zz_0_xz, g_x_0_0_0_yy_zz_0_yy, g_x_0_0_0_yy_zz_0_yz, g_x_0_0_0_yy_zz_0_zz, g_xyy_zz_0_xx, g_xyy_zz_0_xy, g_xyy_zz_0_xz, g_xyy_zz_0_yy, g_xyy_zz_0_yz, g_xyy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_0_xx[i] = 2.0 * g_xyy_zz_0_xx[i] * a_exp;

        g_x_0_0_0_yy_zz_0_xy[i] = 2.0 * g_xyy_zz_0_xy[i] * a_exp;

        g_x_0_0_0_yy_zz_0_xz[i] = 2.0 * g_xyy_zz_0_xz[i] * a_exp;

        g_x_0_0_0_yy_zz_0_yy[i] = 2.0 * g_xyy_zz_0_yy[i] * a_exp;

        g_x_0_0_0_yy_zz_0_yz[i] = 2.0 * g_xyy_zz_0_yz[i] * a_exp;

        g_x_0_0_0_yy_zz_0_zz[i] = 2.0 * g_xyy_zz_0_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_0_xx, g_x_0_0_0_yz_xx_0_xy, g_x_0_0_0_yz_xx_0_xz, g_x_0_0_0_yz_xx_0_yy, g_x_0_0_0_yz_xx_0_yz, g_x_0_0_0_yz_xx_0_zz, g_xyz_xx_0_xx, g_xyz_xx_0_xy, g_xyz_xx_0_xz, g_xyz_xx_0_yy, g_xyz_xx_0_yz, g_xyz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_0_xx[i] = 2.0 * g_xyz_xx_0_xx[i] * a_exp;

        g_x_0_0_0_yz_xx_0_xy[i] = 2.0 * g_xyz_xx_0_xy[i] * a_exp;

        g_x_0_0_0_yz_xx_0_xz[i] = 2.0 * g_xyz_xx_0_xz[i] * a_exp;

        g_x_0_0_0_yz_xx_0_yy[i] = 2.0 * g_xyz_xx_0_yy[i] * a_exp;

        g_x_0_0_0_yz_xx_0_yz[i] = 2.0 * g_xyz_xx_0_yz[i] * a_exp;

        g_x_0_0_0_yz_xx_0_zz[i] = 2.0 * g_xyz_xx_0_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_0_xx, g_x_0_0_0_yz_xy_0_xy, g_x_0_0_0_yz_xy_0_xz, g_x_0_0_0_yz_xy_0_yy, g_x_0_0_0_yz_xy_0_yz, g_x_0_0_0_yz_xy_0_zz, g_xyz_xy_0_xx, g_xyz_xy_0_xy, g_xyz_xy_0_xz, g_xyz_xy_0_yy, g_xyz_xy_0_yz, g_xyz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_0_xx[i] = 2.0 * g_xyz_xy_0_xx[i] * a_exp;

        g_x_0_0_0_yz_xy_0_xy[i] = 2.0 * g_xyz_xy_0_xy[i] * a_exp;

        g_x_0_0_0_yz_xy_0_xz[i] = 2.0 * g_xyz_xy_0_xz[i] * a_exp;

        g_x_0_0_0_yz_xy_0_yy[i] = 2.0 * g_xyz_xy_0_yy[i] * a_exp;

        g_x_0_0_0_yz_xy_0_yz[i] = 2.0 * g_xyz_xy_0_yz[i] * a_exp;

        g_x_0_0_0_yz_xy_0_zz[i] = 2.0 * g_xyz_xy_0_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_0_xx, g_x_0_0_0_yz_xz_0_xy, g_x_0_0_0_yz_xz_0_xz, g_x_0_0_0_yz_xz_0_yy, g_x_0_0_0_yz_xz_0_yz, g_x_0_0_0_yz_xz_0_zz, g_xyz_xz_0_xx, g_xyz_xz_0_xy, g_xyz_xz_0_xz, g_xyz_xz_0_yy, g_xyz_xz_0_yz, g_xyz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_0_xx[i] = 2.0 * g_xyz_xz_0_xx[i] * a_exp;

        g_x_0_0_0_yz_xz_0_xy[i] = 2.0 * g_xyz_xz_0_xy[i] * a_exp;

        g_x_0_0_0_yz_xz_0_xz[i] = 2.0 * g_xyz_xz_0_xz[i] * a_exp;

        g_x_0_0_0_yz_xz_0_yy[i] = 2.0 * g_xyz_xz_0_yy[i] * a_exp;

        g_x_0_0_0_yz_xz_0_yz[i] = 2.0 * g_xyz_xz_0_yz[i] * a_exp;

        g_x_0_0_0_yz_xz_0_zz[i] = 2.0 * g_xyz_xz_0_zz[i] * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_0_xx, g_x_0_0_0_yz_yy_0_xy, g_x_0_0_0_yz_yy_0_xz, g_x_0_0_0_yz_yy_0_yy, g_x_0_0_0_yz_yy_0_yz, g_x_0_0_0_yz_yy_0_zz, g_xyz_yy_0_xx, g_xyz_yy_0_xy, g_xyz_yy_0_xz, g_xyz_yy_0_yy, g_xyz_yy_0_yz, g_xyz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_0_xx[i] = 2.0 * g_xyz_yy_0_xx[i] * a_exp;

        g_x_0_0_0_yz_yy_0_xy[i] = 2.0 * g_xyz_yy_0_xy[i] * a_exp;

        g_x_0_0_0_yz_yy_0_xz[i] = 2.0 * g_xyz_yy_0_xz[i] * a_exp;

        g_x_0_0_0_yz_yy_0_yy[i] = 2.0 * g_xyz_yy_0_yy[i] * a_exp;

        g_x_0_0_0_yz_yy_0_yz[i] = 2.0 * g_xyz_yy_0_yz[i] * a_exp;

        g_x_0_0_0_yz_yy_0_zz[i] = 2.0 * g_xyz_yy_0_zz[i] * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_0_xx, g_x_0_0_0_yz_yz_0_xy, g_x_0_0_0_yz_yz_0_xz, g_x_0_0_0_yz_yz_0_yy, g_x_0_0_0_yz_yz_0_yz, g_x_0_0_0_yz_yz_0_zz, g_xyz_yz_0_xx, g_xyz_yz_0_xy, g_xyz_yz_0_xz, g_xyz_yz_0_yy, g_xyz_yz_0_yz, g_xyz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_0_xx[i] = 2.0 * g_xyz_yz_0_xx[i] * a_exp;

        g_x_0_0_0_yz_yz_0_xy[i] = 2.0 * g_xyz_yz_0_xy[i] * a_exp;

        g_x_0_0_0_yz_yz_0_xz[i] = 2.0 * g_xyz_yz_0_xz[i] * a_exp;

        g_x_0_0_0_yz_yz_0_yy[i] = 2.0 * g_xyz_yz_0_yy[i] * a_exp;

        g_x_0_0_0_yz_yz_0_yz[i] = 2.0 * g_xyz_yz_0_yz[i] * a_exp;

        g_x_0_0_0_yz_yz_0_zz[i] = 2.0 * g_xyz_yz_0_zz[i] * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_0_xx, g_x_0_0_0_yz_zz_0_xy, g_x_0_0_0_yz_zz_0_xz, g_x_0_0_0_yz_zz_0_yy, g_x_0_0_0_yz_zz_0_yz, g_x_0_0_0_yz_zz_0_zz, g_xyz_zz_0_xx, g_xyz_zz_0_xy, g_xyz_zz_0_xz, g_xyz_zz_0_yy, g_xyz_zz_0_yz, g_xyz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_0_xx[i] = 2.0 * g_xyz_zz_0_xx[i] * a_exp;

        g_x_0_0_0_yz_zz_0_xy[i] = 2.0 * g_xyz_zz_0_xy[i] * a_exp;

        g_x_0_0_0_yz_zz_0_xz[i] = 2.0 * g_xyz_zz_0_xz[i] * a_exp;

        g_x_0_0_0_yz_zz_0_yy[i] = 2.0 * g_xyz_zz_0_yy[i] * a_exp;

        g_x_0_0_0_yz_zz_0_yz[i] = 2.0 * g_xyz_zz_0_yz[i] * a_exp;

        g_x_0_0_0_yz_zz_0_zz[i] = 2.0 * g_xyz_zz_0_zz[i] * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_0_xx, g_x_0_0_0_zz_xx_0_xy, g_x_0_0_0_zz_xx_0_xz, g_x_0_0_0_zz_xx_0_yy, g_x_0_0_0_zz_xx_0_yz, g_x_0_0_0_zz_xx_0_zz, g_xzz_xx_0_xx, g_xzz_xx_0_xy, g_xzz_xx_0_xz, g_xzz_xx_0_yy, g_xzz_xx_0_yz, g_xzz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_0_xx[i] = 2.0 * g_xzz_xx_0_xx[i] * a_exp;

        g_x_0_0_0_zz_xx_0_xy[i] = 2.0 * g_xzz_xx_0_xy[i] * a_exp;

        g_x_0_0_0_zz_xx_0_xz[i] = 2.0 * g_xzz_xx_0_xz[i] * a_exp;

        g_x_0_0_0_zz_xx_0_yy[i] = 2.0 * g_xzz_xx_0_yy[i] * a_exp;

        g_x_0_0_0_zz_xx_0_yz[i] = 2.0 * g_xzz_xx_0_yz[i] * a_exp;

        g_x_0_0_0_zz_xx_0_zz[i] = 2.0 * g_xzz_xx_0_zz[i] * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_0_xx, g_x_0_0_0_zz_xy_0_xy, g_x_0_0_0_zz_xy_0_xz, g_x_0_0_0_zz_xy_0_yy, g_x_0_0_0_zz_xy_0_yz, g_x_0_0_0_zz_xy_0_zz, g_xzz_xy_0_xx, g_xzz_xy_0_xy, g_xzz_xy_0_xz, g_xzz_xy_0_yy, g_xzz_xy_0_yz, g_xzz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_0_xx[i] = 2.0 * g_xzz_xy_0_xx[i] * a_exp;

        g_x_0_0_0_zz_xy_0_xy[i] = 2.0 * g_xzz_xy_0_xy[i] * a_exp;

        g_x_0_0_0_zz_xy_0_xz[i] = 2.0 * g_xzz_xy_0_xz[i] * a_exp;

        g_x_0_0_0_zz_xy_0_yy[i] = 2.0 * g_xzz_xy_0_yy[i] * a_exp;

        g_x_0_0_0_zz_xy_0_yz[i] = 2.0 * g_xzz_xy_0_yz[i] * a_exp;

        g_x_0_0_0_zz_xy_0_zz[i] = 2.0 * g_xzz_xy_0_zz[i] * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_0_xx, g_x_0_0_0_zz_xz_0_xy, g_x_0_0_0_zz_xz_0_xz, g_x_0_0_0_zz_xz_0_yy, g_x_0_0_0_zz_xz_0_yz, g_x_0_0_0_zz_xz_0_zz, g_xzz_xz_0_xx, g_xzz_xz_0_xy, g_xzz_xz_0_xz, g_xzz_xz_0_yy, g_xzz_xz_0_yz, g_xzz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_0_xx[i] = 2.0 * g_xzz_xz_0_xx[i] * a_exp;

        g_x_0_0_0_zz_xz_0_xy[i] = 2.0 * g_xzz_xz_0_xy[i] * a_exp;

        g_x_0_0_0_zz_xz_0_xz[i] = 2.0 * g_xzz_xz_0_xz[i] * a_exp;

        g_x_0_0_0_zz_xz_0_yy[i] = 2.0 * g_xzz_xz_0_yy[i] * a_exp;

        g_x_0_0_0_zz_xz_0_yz[i] = 2.0 * g_xzz_xz_0_yz[i] * a_exp;

        g_x_0_0_0_zz_xz_0_zz[i] = 2.0 * g_xzz_xz_0_zz[i] * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_0_xx, g_x_0_0_0_zz_yy_0_xy, g_x_0_0_0_zz_yy_0_xz, g_x_0_0_0_zz_yy_0_yy, g_x_0_0_0_zz_yy_0_yz, g_x_0_0_0_zz_yy_0_zz, g_xzz_yy_0_xx, g_xzz_yy_0_xy, g_xzz_yy_0_xz, g_xzz_yy_0_yy, g_xzz_yy_0_yz, g_xzz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_0_xx[i] = 2.0 * g_xzz_yy_0_xx[i] * a_exp;

        g_x_0_0_0_zz_yy_0_xy[i] = 2.0 * g_xzz_yy_0_xy[i] * a_exp;

        g_x_0_0_0_zz_yy_0_xz[i] = 2.0 * g_xzz_yy_0_xz[i] * a_exp;

        g_x_0_0_0_zz_yy_0_yy[i] = 2.0 * g_xzz_yy_0_yy[i] * a_exp;

        g_x_0_0_0_zz_yy_0_yz[i] = 2.0 * g_xzz_yy_0_yz[i] * a_exp;

        g_x_0_0_0_zz_yy_0_zz[i] = 2.0 * g_xzz_yy_0_zz[i] * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_0_xx, g_x_0_0_0_zz_yz_0_xy, g_x_0_0_0_zz_yz_0_xz, g_x_0_0_0_zz_yz_0_yy, g_x_0_0_0_zz_yz_0_yz, g_x_0_0_0_zz_yz_0_zz, g_xzz_yz_0_xx, g_xzz_yz_0_xy, g_xzz_yz_0_xz, g_xzz_yz_0_yy, g_xzz_yz_0_yz, g_xzz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_0_xx[i] = 2.0 * g_xzz_yz_0_xx[i] * a_exp;

        g_x_0_0_0_zz_yz_0_xy[i] = 2.0 * g_xzz_yz_0_xy[i] * a_exp;

        g_x_0_0_0_zz_yz_0_xz[i] = 2.0 * g_xzz_yz_0_xz[i] * a_exp;

        g_x_0_0_0_zz_yz_0_yy[i] = 2.0 * g_xzz_yz_0_yy[i] * a_exp;

        g_x_0_0_0_zz_yz_0_yz[i] = 2.0 * g_xzz_yz_0_yz[i] * a_exp;

        g_x_0_0_0_zz_yz_0_zz[i] = 2.0 * g_xzz_yz_0_zz[i] * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_0_xx, g_x_0_0_0_zz_zz_0_xy, g_x_0_0_0_zz_zz_0_xz, g_x_0_0_0_zz_zz_0_yy, g_x_0_0_0_zz_zz_0_yz, g_x_0_0_0_zz_zz_0_zz, g_xzz_zz_0_xx, g_xzz_zz_0_xy, g_xzz_zz_0_xz, g_xzz_zz_0_yy, g_xzz_zz_0_yz, g_xzz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_0_xx[i] = 2.0 * g_xzz_zz_0_xx[i] * a_exp;

        g_x_0_0_0_zz_zz_0_xy[i] = 2.0 * g_xzz_zz_0_xy[i] * a_exp;

        g_x_0_0_0_zz_zz_0_xz[i] = 2.0 * g_xzz_zz_0_xz[i] * a_exp;

        g_x_0_0_0_zz_zz_0_yy[i] = 2.0 * g_xzz_zz_0_yy[i] * a_exp;

        g_x_0_0_0_zz_zz_0_yz[i] = 2.0 * g_xzz_zz_0_yz[i] * a_exp;

        g_x_0_0_0_zz_zz_0_zz[i] = 2.0 * g_xzz_zz_0_zz[i] * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xxy_xx_0_xx, g_xxy_xx_0_xy, g_xxy_xx_0_xz, g_xxy_xx_0_yy, g_xxy_xx_0_yz, g_xxy_xx_0_zz, g_y_0_0_0_xx_xx_0_xx, g_y_0_0_0_xx_xx_0_xy, g_y_0_0_0_xx_xx_0_xz, g_y_0_0_0_xx_xx_0_yy, g_y_0_0_0_xx_xx_0_yz, g_y_0_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_0_xx[i] = 2.0 * g_xxy_xx_0_xx[i] * a_exp;

        g_y_0_0_0_xx_xx_0_xy[i] = 2.0 * g_xxy_xx_0_xy[i] * a_exp;

        g_y_0_0_0_xx_xx_0_xz[i] = 2.0 * g_xxy_xx_0_xz[i] * a_exp;

        g_y_0_0_0_xx_xx_0_yy[i] = 2.0 * g_xxy_xx_0_yy[i] * a_exp;

        g_y_0_0_0_xx_xx_0_yz[i] = 2.0 * g_xxy_xx_0_yz[i] * a_exp;

        g_y_0_0_0_xx_xx_0_zz[i] = 2.0 * g_xxy_xx_0_zz[i] * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xxy_xy_0_xx, g_xxy_xy_0_xy, g_xxy_xy_0_xz, g_xxy_xy_0_yy, g_xxy_xy_0_yz, g_xxy_xy_0_zz, g_y_0_0_0_xx_xy_0_xx, g_y_0_0_0_xx_xy_0_xy, g_y_0_0_0_xx_xy_0_xz, g_y_0_0_0_xx_xy_0_yy, g_y_0_0_0_xx_xy_0_yz, g_y_0_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_0_xx[i] = 2.0 * g_xxy_xy_0_xx[i] * a_exp;

        g_y_0_0_0_xx_xy_0_xy[i] = 2.0 * g_xxy_xy_0_xy[i] * a_exp;

        g_y_0_0_0_xx_xy_0_xz[i] = 2.0 * g_xxy_xy_0_xz[i] * a_exp;

        g_y_0_0_0_xx_xy_0_yy[i] = 2.0 * g_xxy_xy_0_yy[i] * a_exp;

        g_y_0_0_0_xx_xy_0_yz[i] = 2.0 * g_xxy_xy_0_yz[i] * a_exp;

        g_y_0_0_0_xx_xy_0_zz[i] = 2.0 * g_xxy_xy_0_zz[i] * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xxy_xz_0_xx, g_xxy_xz_0_xy, g_xxy_xz_0_xz, g_xxy_xz_0_yy, g_xxy_xz_0_yz, g_xxy_xz_0_zz, g_y_0_0_0_xx_xz_0_xx, g_y_0_0_0_xx_xz_0_xy, g_y_0_0_0_xx_xz_0_xz, g_y_0_0_0_xx_xz_0_yy, g_y_0_0_0_xx_xz_0_yz, g_y_0_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_0_xx[i] = 2.0 * g_xxy_xz_0_xx[i] * a_exp;

        g_y_0_0_0_xx_xz_0_xy[i] = 2.0 * g_xxy_xz_0_xy[i] * a_exp;

        g_y_0_0_0_xx_xz_0_xz[i] = 2.0 * g_xxy_xz_0_xz[i] * a_exp;

        g_y_0_0_0_xx_xz_0_yy[i] = 2.0 * g_xxy_xz_0_yy[i] * a_exp;

        g_y_0_0_0_xx_xz_0_yz[i] = 2.0 * g_xxy_xz_0_yz[i] * a_exp;

        g_y_0_0_0_xx_xz_0_zz[i] = 2.0 * g_xxy_xz_0_zz[i] * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_xxy_yy_0_xx, g_xxy_yy_0_xy, g_xxy_yy_0_xz, g_xxy_yy_0_yy, g_xxy_yy_0_yz, g_xxy_yy_0_zz, g_y_0_0_0_xx_yy_0_xx, g_y_0_0_0_xx_yy_0_xy, g_y_0_0_0_xx_yy_0_xz, g_y_0_0_0_xx_yy_0_yy, g_y_0_0_0_xx_yy_0_yz, g_y_0_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_0_xx[i] = 2.0 * g_xxy_yy_0_xx[i] * a_exp;

        g_y_0_0_0_xx_yy_0_xy[i] = 2.0 * g_xxy_yy_0_xy[i] * a_exp;

        g_y_0_0_0_xx_yy_0_xz[i] = 2.0 * g_xxy_yy_0_xz[i] * a_exp;

        g_y_0_0_0_xx_yy_0_yy[i] = 2.0 * g_xxy_yy_0_yy[i] * a_exp;

        g_y_0_0_0_xx_yy_0_yz[i] = 2.0 * g_xxy_yy_0_yz[i] * a_exp;

        g_y_0_0_0_xx_yy_0_zz[i] = 2.0 * g_xxy_yy_0_zz[i] * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_xxy_yz_0_xx, g_xxy_yz_0_xy, g_xxy_yz_0_xz, g_xxy_yz_0_yy, g_xxy_yz_0_yz, g_xxy_yz_0_zz, g_y_0_0_0_xx_yz_0_xx, g_y_0_0_0_xx_yz_0_xy, g_y_0_0_0_xx_yz_0_xz, g_y_0_0_0_xx_yz_0_yy, g_y_0_0_0_xx_yz_0_yz, g_y_0_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_0_xx[i] = 2.0 * g_xxy_yz_0_xx[i] * a_exp;

        g_y_0_0_0_xx_yz_0_xy[i] = 2.0 * g_xxy_yz_0_xy[i] * a_exp;

        g_y_0_0_0_xx_yz_0_xz[i] = 2.0 * g_xxy_yz_0_xz[i] * a_exp;

        g_y_0_0_0_xx_yz_0_yy[i] = 2.0 * g_xxy_yz_0_yy[i] * a_exp;

        g_y_0_0_0_xx_yz_0_yz[i] = 2.0 * g_xxy_yz_0_yz[i] * a_exp;

        g_y_0_0_0_xx_yz_0_zz[i] = 2.0 * g_xxy_yz_0_zz[i] * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_xxy_zz_0_xx, g_xxy_zz_0_xy, g_xxy_zz_0_xz, g_xxy_zz_0_yy, g_xxy_zz_0_yz, g_xxy_zz_0_zz, g_y_0_0_0_xx_zz_0_xx, g_y_0_0_0_xx_zz_0_xy, g_y_0_0_0_xx_zz_0_xz, g_y_0_0_0_xx_zz_0_yy, g_y_0_0_0_xx_zz_0_yz, g_y_0_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_0_xx[i] = 2.0 * g_xxy_zz_0_xx[i] * a_exp;

        g_y_0_0_0_xx_zz_0_xy[i] = 2.0 * g_xxy_zz_0_xy[i] * a_exp;

        g_y_0_0_0_xx_zz_0_xz[i] = 2.0 * g_xxy_zz_0_xz[i] * a_exp;

        g_y_0_0_0_xx_zz_0_yy[i] = 2.0 * g_xxy_zz_0_yy[i] * a_exp;

        g_y_0_0_0_xx_zz_0_yz[i] = 2.0 * g_xxy_zz_0_yz[i] * a_exp;

        g_y_0_0_0_xx_zz_0_zz[i] = 2.0 * g_xxy_zz_0_zz[i] * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_xyy_xx_0_xx, g_xyy_xx_0_xy, g_xyy_xx_0_xz, g_xyy_xx_0_yy, g_xyy_xx_0_yz, g_xyy_xx_0_zz, g_y_0_0_0_xy_xx_0_xx, g_y_0_0_0_xy_xx_0_xy, g_y_0_0_0_xy_xx_0_xz, g_y_0_0_0_xy_xx_0_yy, g_y_0_0_0_xy_xx_0_yz, g_y_0_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_0_xx[i] = -g_x_xx_0_xx[i] + 2.0 * g_xyy_xx_0_xx[i] * a_exp;

        g_y_0_0_0_xy_xx_0_xy[i] = -g_x_xx_0_xy[i] + 2.0 * g_xyy_xx_0_xy[i] * a_exp;

        g_y_0_0_0_xy_xx_0_xz[i] = -g_x_xx_0_xz[i] + 2.0 * g_xyy_xx_0_xz[i] * a_exp;

        g_y_0_0_0_xy_xx_0_yy[i] = -g_x_xx_0_yy[i] + 2.0 * g_xyy_xx_0_yy[i] * a_exp;

        g_y_0_0_0_xy_xx_0_yz[i] = -g_x_xx_0_yz[i] + 2.0 * g_xyy_xx_0_yz[i] * a_exp;

        g_y_0_0_0_xy_xx_0_zz[i] = -g_x_xx_0_zz[i] + 2.0 * g_xyy_xx_0_zz[i] * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_xyy_xy_0_xx, g_xyy_xy_0_xy, g_xyy_xy_0_xz, g_xyy_xy_0_yy, g_xyy_xy_0_yz, g_xyy_xy_0_zz, g_y_0_0_0_xy_xy_0_xx, g_y_0_0_0_xy_xy_0_xy, g_y_0_0_0_xy_xy_0_xz, g_y_0_0_0_xy_xy_0_yy, g_y_0_0_0_xy_xy_0_yz, g_y_0_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_0_xx[i] = -g_x_xy_0_xx[i] + 2.0 * g_xyy_xy_0_xx[i] * a_exp;

        g_y_0_0_0_xy_xy_0_xy[i] = -g_x_xy_0_xy[i] + 2.0 * g_xyy_xy_0_xy[i] * a_exp;

        g_y_0_0_0_xy_xy_0_xz[i] = -g_x_xy_0_xz[i] + 2.0 * g_xyy_xy_0_xz[i] * a_exp;

        g_y_0_0_0_xy_xy_0_yy[i] = -g_x_xy_0_yy[i] + 2.0 * g_xyy_xy_0_yy[i] * a_exp;

        g_y_0_0_0_xy_xy_0_yz[i] = -g_x_xy_0_yz[i] + 2.0 * g_xyy_xy_0_yz[i] * a_exp;

        g_y_0_0_0_xy_xy_0_zz[i] = -g_x_xy_0_zz[i] + 2.0 * g_xyy_xy_0_zz[i] * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_xyy_xz_0_xx, g_xyy_xz_0_xy, g_xyy_xz_0_xz, g_xyy_xz_0_yy, g_xyy_xz_0_yz, g_xyy_xz_0_zz, g_y_0_0_0_xy_xz_0_xx, g_y_0_0_0_xy_xz_0_xy, g_y_0_0_0_xy_xz_0_xz, g_y_0_0_0_xy_xz_0_yy, g_y_0_0_0_xy_xz_0_yz, g_y_0_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_0_xx[i] = -g_x_xz_0_xx[i] + 2.0 * g_xyy_xz_0_xx[i] * a_exp;

        g_y_0_0_0_xy_xz_0_xy[i] = -g_x_xz_0_xy[i] + 2.0 * g_xyy_xz_0_xy[i] * a_exp;

        g_y_0_0_0_xy_xz_0_xz[i] = -g_x_xz_0_xz[i] + 2.0 * g_xyy_xz_0_xz[i] * a_exp;

        g_y_0_0_0_xy_xz_0_yy[i] = -g_x_xz_0_yy[i] + 2.0 * g_xyy_xz_0_yy[i] * a_exp;

        g_y_0_0_0_xy_xz_0_yz[i] = -g_x_xz_0_yz[i] + 2.0 * g_xyy_xz_0_yz[i] * a_exp;

        g_y_0_0_0_xy_xz_0_zz[i] = -g_x_xz_0_zz[i] + 2.0 * g_xyy_xz_0_zz[i] * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_xyy_yy_0_xx, g_xyy_yy_0_xy, g_xyy_yy_0_xz, g_xyy_yy_0_yy, g_xyy_yy_0_yz, g_xyy_yy_0_zz, g_y_0_0_0_xy_yy_0_xx, g_y_0_0_0_xy_yy_0_xy, g_y_0_0_0_xy_yy_0_xz, g_y_0_0_0_xy_yy_0_yy, g_y_0_0_0_xy_yy_0_yz, g_y_0_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_0_xx[i] = -g_x_yy_0_xx[i] + 2.0 * g_xyy_yy_0_xx[i] * a_exp;

        g_y_0_0_0_xy_yy_0_xy[i] = -g_x_yy_0_xy[i] + 2.0 * g_xyy_yy_0_xy[i] * a_exp;

        g_y_0_0_0_xy_yy_0_xz[i] = -g_x_yy_0_xz[i] + 2.0 * g_xyy_yy_0_xz[i] * a_exp;

        g_y_0_0_0_xy_yy_0_yy[i] = -g_x_yy_0_yy[i] + 2.0 * g_xyy_yy_0_yy[i] * a_exp;

        g_y_0_0_0_xy_yy_0_yz[i] = -g_x_yy_0_yz[i] + 2.0 * g_xyy_yy_0_yz[i] * a_exp;

        g_y_0_0_0_xy_yy_0_zz[i] = -g_x_yy_0_zz[i] + 2.0 * g_xyy_yy_0_zz[i] * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_xyy_yz_0_xx, g_xyy_yz_0_xy, g_xyy_yz_0_xz, g_xyy_yz_0_yy, g_xyy_yz_0_yz, g_xyy_yz_0_zz, g_y_0_0_0_xy_yz_0_xx, g_y_0_0_0_xy_yz_0_xy, g_y_0_0_0_xy_yz_0_xz, g_y_0_0_0_xy_yz_0_yy, g_y_0_0_0_xy_yz_0_yz, g_y_0_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_0_xx[i] = -g_x_yz_0_xx[i] + 2.0 * g_xyy_yz_0_xx[i] * a_exp;

        g_y_0_0_0_xy_yz_0_xy[i] = -g_x_yz_0_xy[i] + 2.0 * g_xyy_yz_0_xy[i] * a_exp;

        g_y_0_0_0_xy_yz_0_xz[i] = -g_x_yz_0_xz[i] + 2.0 * g_xyy_yz_0_xz[i] * a_exp;

        g_y_0_0_0_xy_yz_0_yy[i] = -g_x_yz_0_yy[i] + 2.0 * g_xyy_yz_0_yy[i] * a_exp;

        g_y_0_0_0_xy_yz_0_yz[i] = -g_x_yz_0_yz[i] + 2.0 * g_xyy_yz_0_yz[i] * a_exp;

        g_y_0_0_0_xy_yz_0_zz[i] = -g_x_yz_0_zz[i] + 2.0 * g_xyy_yz_0_zz[i] * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_xyy_zz_0_xx, g_xyy_zz_0_xy, g_xyy_zz_0_xz, g_xyy_zz_0_yy, g_xyy_zz_0_yz, g_xyy_zz_0_zz, g_y_0_0_0_xy_zz_0_xx, g_y_0_0_0_xy_zz_0_xy, g_y_0_0_0_xy_zz_0_xz, g_y_0_0_0_xy_zz_0_yy, g_y_0_0_0_xy_zz_0_yz, g_y_0_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_0_xx[i] = -g_x_zz_0_xx[i] + 2.0 * g_xyy_zz_0_xx[i] * a_exp;

        g_y_0_0_0_xy_zz_0_xy[i] = -g_x_zz_0_xy[i] + 2.0 * g_xyy_zz_0_xy[i] * a_exp;

        g_y_0_0_0_xy_zz_0_xz[i] = -g_x_zz_0_xz[i] + 2.0 * g_xyy_zz_0_xz[i] * a_exp;

        g_y_0_0_0_xy_zz_0_yy[i] = -g_x_zz_0_yy[i] + 2.0 * g_xyy_zz_0_yy[i] * a_exp;

        g_y_0_0_0_xy_zz_0_yz[i] = -g_x_zz_0_yz[i] + 2.0 * g_xyy_zz_0_yz[i] * a_exp;

        g_y_0_0_0_xy_zz_0_zz[i] = -g_x_zz_0_zz[i] + 2.0 * g_xyy_zz_0_zz[i] * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_xyz_xx_0_xx, g_xyz_xx_0_xy, g_xyz_xx_0_xz, g_xyz_xx_0_yy, g_xyz_xx_0_yz, g_xyz_xx_0_zz, g_y_0_0_0_xz_xx_0_xx, g_y_0_0_0_xz_xx_0_xy, g_y_0_0_0_xz_xx_0_xz, g_y_0_0_0_xz_xx_0_yy, g_y_0_0_0_xz_xx_0_yz, g_y_0_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_0_xx[i] = 2.0 * g_xyz_xx_0_xx[i] * a_exp;

        g_y_0_0_0_xz_xx_0_xy[i] = 2.0 * g_xyz_xx_0_xy[i] * a_exp;

        g_y_0_0_0_xz_xx_0_xz[i] = 2.0 * g_xyz_xx_0_xz[i] * a_exp;

        g_y_0_0_0_xz_xx_0_yy[i] = 2.0 * g_xyz_xx_0_yy[i] * a_exp;

        g_y_0_0_0_xz_xx_0_yz[i] = 2.0 * g_xyz_xx_0_yz[i] * a_exp;

        g_y_0_0_0_xz_xx_0_zz[i] = 2.0 * g_xyz_xx_0_zz[i] * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_xyz_xy_0_xx, g_xyz_xy_0_xy, g_xyz_xy_0_xz, g_xyz_xy_0_yy, g_xyz_xy_0_yz, g_xyz_xy_0_zz, g_y_0_0_0_xz_xy_0_xx, g_y_0_0_0_xz_xy_0_xy, g_y_0_0_0_xz_xy_0_xz, g_y_0_0_0_xz_xy_0_yy, g_y_0_0_0_xz_xy_0_yz, g_y_0_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_0_xx[i] = 2.0 * g_xyz_xy_0_xx[i] * a_exp;

        g_y_0_0_0_xz_xy_0_xy[i] = 2.0 * g_xyz_xy_0_xy[i] * a_exp;

        g_y_0_0_0_xz_xy_0_xz[i] = 2.0 * g_xyz_xy_0_xz[i] * a_exp;

        g_y_0_0_0_xz_xy_0_yy[i] = 2.0 * g_xyz_xy_0_yy[i] * a_exp;

        g_y_0_0_0_xz_xy_0_yz[i] = 2.0 * g_xyz_xy_0_yz[i] * a_exp;

        g_y_0_0_0_xz_xy_0_zz[i] = 2.0 * g_xyz_xy_0_zz[i] * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_xyz_xz_0_xx, g_xyz_xz_0_xy, g_xyz_xz_0_xz, g_xyz_xz_0_yy, g_xyz_xz_0_yz, g_xyz_xz_0_zz, g_y_0_0_0_xz_xz_0_xx, g_y_0_0_0_xz_xz_0_xy, g_y_0_0_0_xz_xz_0_xz, g_y_0_0_0_xz_xz_0_yy, g_y_0_0_0_xz_xz_0_yz, g_y_0_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_0_xx[i] = 2.0 * g_xyz_xz_0_xx[i] * a_exp;

        g_y_0_0_0_xz_xz_0_xy[i] = 2.0 * g_xyz_xz_0_xy[i] * a_exp;

        g_y_0_0_0_xz_xz_0_xz[i] = 2.0 * g_xyz_xz_0_xz[i] * a_exp;

        g_y_0_0_0_xz_xz_0_yy[i] = 2.0 * g_xyz_xz_0_yy[i] * a_exp;

        g_y_0_0_0_xz_xz_0_yz[i] = 2.0 * g_xyz_xz_0_yz[i] * a_exp;

        g_y_0_0_0_xz_xz_0_zz[i] = 2.0 * g_xyz_xz_0_zz[i] * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_xyz_yy_0_xx, g_xyz_yy_0_xy, g_xyz_yy_0_xz, g_xyz_yy_0_yy, g_xyz_yy_0_yz, g_xyz_yy_0_zz, g_y_0_0_0_xz_yy_0_xx, g_y_0_0_0_xz_yy_0_xy, g_y_0_0_0_xz_yy_0_xz, g_y_0_0_0_xz_yy_0_yy, g_y_0_0_0_xz_yy_0_yz, g_y_0_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_0_xx[i] = 2.0 * g_xyz_yy_0_xx[i] * a_exp;

        g_y_0_0_0_xz_yy_0_xy[i] = 2.0 * g_xyz_yy_0_xy[i] * a_exp;

        g_y_0_0_0_xz_yy_0_xz[i] = 2.0 * g_xyz_yy_0_xz[i] * a_exp;

        g_y_0_0_0_xz_yy_0_yy[i] = 2.0 * g_xyz_yy_0_yy[i] * a_exp;

        g_y_0_0_0_xz_yy_0_yz[i] = 2.0 * g_xyz_yy_0_yz[i] * a_exp;

        g_y_0_0_0_xz_yy_0_zz[i] = 2.0 * g_xyz_yy_0_zz[i] * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_xyz_yz_0_xx, g_xyz_yz_0_xy, g_xyz_yz_0_xz, g_xyz_yz_0_yy, g_xyz_yz_0_yz, g_xyz_yz_0_zz, g_y_0_0_0_xz_yz_0_xx, g_y_0_0_0_xz_yz_0_xy, g_y_0_0_0_xz_yz_0_xz, g_y_0_0_0_xz_yz_0_yy, g_y_0_0_0_xz_yz_0_yz, g_y_0_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_0_xx[i] = 2.0 * g_xyz_yz_0_xx[i] * a_exp;

        g_y_0_0_0_xz_yz_0_xy[i] = 2.0 * g_xyz_yz_0_xy[i] * a_exp;

        g_y_0_0_0_xz_yz_0_xz[i] = 2.0 * g_xyz_yz_0_xz[i] * a_exp;

        g_y_0_0_0_xz_yz_0_yy[i] = 2.0 * g_xyz_yz_0_yy[i] * a_exp;

        g_y_0_0_0_xz_yz_0_yz[i] = 2.0 * g_xyz_yz_0_yz[i] * a_exp;

        g_y_0_0_0_xz_yz_0_zz[i] = 2.0 * g_xyz_yz_0_zz[i] * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_xyz_zz_0_xx, g_xyz_zz_0_xy, g_xyz_zz_0_xz, g_xyz_zz_0_yy, g_xyz_zz_0_yz, g_xyz_zz_0_zz, g_y_0_0_0_xz_zz_0_xx, g_y_0_0_0_xz_zz_0_xy, g_y_0_0_0_xz_zz_0_xz, g_y_0_0_0_xz_zz_0_yy, g_y_0_0_0_xz_zz_0_yz, g_y_0_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_0_xx[i] = 2.0 * g_xyz_zz_0_xx[i] * a_exp;

        g_y_0_0_0_xz_zz_0_xy[i] = 2.0 * g_xyz_zz_0_xy[i] * a_exp;

        g_y_0_0_0_xz_zz_0_xz[i] = 2.0 * g_xyz_zz_0_xz[i] * a_exp;

        g_y_0_0_0_xz_zz_0_yy[i] = 2.0 * g_xyz_zz_0_yy[i] * a_exp;

        g_y_0_0_0_xz_zz_0_yz[i] = 2.0 * g_xyz_zz_0_yz[i] * a_exp;

        g_y_0_0_0_xz_zz_0_zz[i] = 2.0 * g_xyz_zz_0_zz[i] * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_0_xx, g_y_0_0_0_yy_xx_0_xy, g_y_0_0_0_yy_xx_0_xz, g_y_0_0_0_yy_xx_0_yy, g_y_0_0_0_yy_xx_0_yz, g_y_0_0_0_yy_xx_0_zz, g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz, g_yyy_xx_0_xx, g_yyy_xx_0_xy, g_yyy_xx_0_xz, g_yyy_xx_0_yy, g_yyy_xx_0_yz, g_yyy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_0_xx[i] = -2.0 * g_y_xx_0_xx[i] + 2.0 * g_yyy_xx_0_xx[i] * a_exp;

        g_y_0_0_0_yy_xx_0_xy[i] = -2.0 * g_y_xx_0_xy[i] + 2.0 * g_yyy_xx_0_xy[i] * a_exp;

        g_y_0_0_0_yy_xx_0_xz[i] = -2.0 * g_y_xx_0_xz[i] + 2.0 * g_yyy_xx_0_xz[i] * a_exp;

        g_y_0_0_0_yy_xx_0_yy[i] = -2.0 * g_y_xx_0_yy[i] + 2.0 * g_yyy_xx_0_yy[i] * a_exp;

        g_y_0_0_0_yy_xx_0_yz[i] = -2.0 * g_y_xx_0_yz[i] + 2.0 * g_yyy_xx_0_yz[i] * a_exp;

        g_y_0_0_0_yy_xx_0_zz[i] = -2.0 * g_y_xx_0_zz[i] + 2.0 * g_yyy_xx_0_zz[i] * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_0_xx, g_y_0_0_0_yy_xy_0_xy, g_y_0_0_0_yy_xy_0_xz, g_y_0_0_0_yy_xy_0_yy, g_y_0_0_0_yy_xy_0_yz, g_y_0_0_0_yy_xy_0_zz, g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz, g_yyy_xy_0_xx, g_yyy_xy_0_xy, g_yyy_xy_0_xz, g_yyy_xy_0_yy, g_yyy_xy_0_yz, g_yyy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_0_xx[i] = -2.0 * g_y_xy_0_xx[i] + 2.0 * g_yyy_xy_0_xx[i] * a_exp;

        g_y_0_0_0_yy_xy_0_xy[i] = -2.0 * g_y_xy_0_xy[i] + 2.0 * g_yyy_xy_0_xy[i] * a_exp;

        g_y_0_0_0_yy_xy_0_xz[i] = -2.0 * g_y_xy_0_xz[i] + 2.0 * g_yyy_xy_0_xz[i] * a_exp;

        g_y_0_0_0_yy_xy_0_yy[i] = -2.0 * g_y_xy_0_yy[i] + 2.0 * g_yyy_xy_0_yy[i] * a_exp;

        g_y_0_0_0_yy_xy_0_yz[i] = -2.0 * g_y_xy_0_yz[i] + 2.0 * g_yyy_xy_0_yz[i] * a_exp;

        g_y_0_0_0_yy_xy_0_zz[i] = -2.0 * g_y_xy_0_zz[i] + 2.0 * g_yyy_xy_0_zz[i] * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_0_xx, g_y_0_0_0_yy_xz_0_xy, g_y_0_0_0_yy_xz_0_xz, g_y_0_0_0_yy_xz_0_yy, g_y_0_0_0_yy_xz_0_yz, g_y_0_0_0_yy_xz_0_zz, g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz, g_yyy_xz_0_xx, g_yyy_xz_0_xy, g_yyy_xz_0_xz, g_yyy_xz_0_yy, g_yyy_xz_0_yz, g_yyy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_0_xx[i] = -2.0 * g_y_xz_0_xx[i] + 2.0 * g_yyy_xz_0_xx[i] * a_exp;

        g_y_0_0_0_yy_xz_0_xy[i] = -2.0 * g_y_xz_0_xy[i] + 2.0 * g_yyy_xz_0_xy[i] * a_exp;

        g_y_0_0_0_yy_xz_0_xz[i] = -2.0 * g_y_xz_0_xz[i] + 2.0 * g_yyy_xz_0_xz[i] * a_exp;

        g_y_0_0_0_yy_xz_0_yy[i] = -2.0 * g_y_xz_0_yy[i] + 2.0 * g_yyy_xz_0_yy[i] * a_exp;

        g_y_0_0_0_yy_xz_0_yz[i] = -2.0 * g_y_xz_0_yz[i] + 2.0 * g_yyy_xz_0_yz[i] * a_exp;

        g_y_0_0_0_yy_xz_0_zz[i] = -2.0 * g_y_xz_0_zz[i] + 2.0 * g_yyy_xz_0_zz[i] * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_0_xx, g_y_0_0_0_yy_yy_0_xy, g_y_0_0_0_yy_yy_0_xz, g_y_0_0_0_yy_yy_0_yy, g_y_0_0_0_yy_yy_0_yz, g_y_0_0_0_yy_yy_0_zz, g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz, g_yyy_yy_0_xx, g_yyy_yy_0_xy, g_yyy_yy_0_xz, g_yyy_yy_0_yy, g_yyy_yy_0_yz, g_yyy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_0_xx[i] = -2.0 * g_y_yy_0_xx[i] + 2.0 * g_yyy_yy_0_xx[i] * a_exp;

        g_y_0_0_0_yy_yy_0_xy[i] = -2.0 * g_y_yy_0_xy[i] + 2.0 * g_yyy_yy_0_xy[i] * a_exp;

        g_y_0_0_0_yy_yy_0_xz[i] = -2.0 * g_y_yy_0_xz[i] + 2.0 * g_yyy_yy_0_xz[i] * a_exp;

        g_y_0_0_0_yy_yy_0_yy[i] = -2.0 * g_y_yy_0_yy[i] + 2.0 * g_yyy_yy_0_yy[i] * a_exp;

        g_y_0_0_0_yy_yy_0_yz[i] = -2.0 * g_y_yy_0_yz[i] + 2.0 * g_yyy_yy_0_yz[i] * a_exp;

        g_y_0_0_0_yy_yy_0_zz[i] = -2.0 * g_y_yy_0_zz[i] + 2.0 * g_yyy_yy_0_zz[i] * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_0_xx, g_y_0_0_0_yy_yz_0_xy, g_y_0_0_0_yy_yz_0_xz, g_y_0_0_0_yy_yz_0_yy, g_y_0_0_0_yy_yz_0_yz, g_y_0_0_0_yy_yz_0_zz, g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz, g_yyy_yz_0_xx, g_yyy_yz_0_xy, g_yyy_yz_0_xz, g_yyy_yz_0_yy, g_yyy_yz_0_yz, g_yyy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_0_xx[i] = -2.0 * g_y_yz_0_xx[i] + 2.0 * g_yyy_yz_0_xx[i] * a_exp;

        g_y_0_0_0_yy_yz_0_xy[i] = -2.0 * g_y_yz_0_xy[i] + 2.0 * g_yyy_yz_0_xy[i] * a_exp;

        g_y_0_0_0_yy_yz_0_xz[i] = -2.0 * g_y_yz_0_xz[i] + 2.0 * g_yyy_yz_0_xz[i] * a_exp;

        g_y_0_0_0_yy_yz_0_yy[i] = -2.0 * g_y_yz_0_yy[i] + 2.0 * g_yyy_yz_0_yy[i] * a_exp;

        g_y_0_0_0_yy_yz_0_yz[i] = -2.0 * g_y_yz_0_yz[i] + 2.0 * g_yyy_yz_0_yz[i] * a_exp;

        g_y_0_0_0_yy_yz_0_zz[i] = -2.0 * g_y_yz_0_zz[i] + 2.0 * g_yyy_yz_0_zz[i] * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_0_xx, g_y_0_0_0_yy_zz_0_xy, g_y_0_0_0_yy_zz_0_xz, g_y_0_0_0_yy_zz_0_yy, g_y_0_0_0_yy_zz_0_yz, g_y_0_0_0_yy_zz_0_zz, g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz, g_yyy_zz_0_xx, g_yyy_zz_0_xy, g_yyy_zz_0_xz, g_yyy_zz_0_yy, g_yyy_zz_0_yz, g_yyy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_0_xx[i] = -2.0 * g_y_zz_0_xx[i] + 2.0 * g_yyy_zz_0_xx[i] * a_exp;

        g_y_0_0_0_yy_zz_0_xy[i] = -2.0 * g_y_zz_0_xy[i] + 2.0 * g_yyy_zz_0_xy[i] * a_exp;

        g_y_0_0_0_yy_zz_0_xz[i] = -2.0 * g_y_zz_0_xz[i] + 2.0 * g_yyy_zz_0_xz[i] * a_exp;

        g_y_0_0_0_yy_zz_0_yy[i] = -2.0 * g_y_zz_0_yy[i] + 2.0 * g_yyy_zz_0_yy[i] * a_exp;

        g_y_0_0_0_yy_zz_0_yz[i] = -2.0 * g_y_zz_0_yz[i] + 2.0 * g_yyy_zz_0_yz[i] * a_exp;

        g_y_0_0_0_yy_zz_0_zz[i] = -2.0 * g_y_zz_0_zz[i] + 2.0 * g_yyy_zz_0_zz[i] * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_0_xx, g_y_0_0_0_yz_xx_0_xy, g_y_0_0_0_yz_xx_0_xz, g_y_0_0_0_yz_xx_0_yy, g_y_0_0_0_yz_xx_0_yz, g_y_0_0_0_yz_xx_0_zz, g_yyz_xx_0_xx, g_yyz_xx_0_xy, g_yyz_xx_0_xz, g_yyz_xx_0_yy, g_yyz_xx_0_yz, g_yyz_xx_0_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_0_xx[i] = -g_z_xx_0_xx[i] + 2.0 * g_yyz_xx_0_xx[i] * a_exp;

        g_y_0_0_0_yz_xx_0_xy[i] = -g_z_xx_0_xy[i] + 2.0 * g_yyz_xx_0_xy[i] * a_exp;

        g_y_0_0_0_yz_xx_0_xz[i] = -g_z_xx_0_xz[i] + 2.0 * g_yyz_xx_0_xz[i] * a_exp;

        g_y_0_0_0_yz_xx_0_yy[i] = -g_z_xx_0_yy[i] + 2.0 * g_yyz_xx_0_yy[i] * a_exp;

        g_y_0_0_0_yz_xx_0_yz[i] = -g_z_xx_0_yz[i] + 2.0 * g_yyz_xx_0_yz[i] * a_exp;

        g_y_0_0_0_yz_xx_0_zz[i] = -g_z_xx_0_zz[i] + 2.0 * g_yyz_xx_0_zz[i] * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_0_xx, g_y_0_0_0_yz_xy_0_xy, g_y_0_0_0_yz_xy_0_xz, g_y_0_0_0_yz_xy_0_yy, g_y_0_0_0_yz_xy_0_yz, g_y_0_0_0_yz_xy_0_zz, g_yyz_xy_0_xx, g_yyz_xy_0_xy, g_yyz_xy_0_xz, g_yyz_xy_0_yy, g_yyz_xy_0_yz, g_yyz_xy_0_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_0_xx[i] = -g_z_xy_0_xx[i] + 2.0 * g_yyz_xy_0_xx[i] * a_exp;

        g_y_0_0_0_yz_xy_0_xy[i] = -g_z_xy_0_xy[i] + 2.0 * g_yyz_xy_0_xy[i] * a_exp;

        g_y_0_0_0_yz_xy_0_xz[i] = -g_z_xy_0_xz[i] + 2.0 * g_yyz_xy_0_xz[i] * a_exp;

        g_y_0_0_0_yz_xy_0_yy[i] = -g_z_xy_0_yy[i] + 2.0 * g_yyz_xy_0_yy[i] * a_exp;

        g_y_0_0_0_yz_xy_0_yz[i] = -g_z_xy_0_yz[i] + 2.0 * g_yyz_xy_0_yz[i] * a_exp;

        g_y_0_0_0_yz_xy_0_zz[i] = -g_z_xy_0_zz[i] + 2.0 * g_yyz_xy_0_zz[i] * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_0_xx, g_y_0_0_0_yz_xz_0_xy, g_y_0_0_0_yz_xz_0_xz, g_y_0_0_0_yz_xz_0_yy, g_y_0_0_0_yz_xz_0_yz, g_y_0_0_0_yz_xz_0_zz, g_yyz_xz_0_xx, g_yyz_xz_0_xy, g_yyz_xz_0_xz, g_yyz_xz_0_yy, g_yyz_xz_0_yz, g_yyz_xz_0_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_0_xx[i] = -g_z_xz_0_xx[i] + 2.0 * g_yyz_xz_0_xx[i] * a_exp;

        g_y_0_0_0_yz_xz_0_xy[i] = -g_z_xz_0_xy[i] + 2.0 * g_yyz_xz_0_xy[i] * a_exp;

        g_y_0_0_0_yz_xz_0_xz[i] = -g_z_xz_0_xz[i] + 2.0 * g_yyz_xz_0_xz[i] * a_exp;

        g_y_0_0_0_yz_xz_0_yy[i] = -g_z_xz_0_yy[i] + 2.0 * g_yyz_xz_0_yy[i] * a_exp;

        g_y_0_0_0_yz_xz_0_yz[i] = -g_z_xz_0_yz[i] + 2.0 * g_yyz_xz_0_yz[i] * a_exp;

        g_y_0_0_0_yz_xz_0_zz[i] = -g_z_xz_0_zz[i] + 2.0 * g_yyz_xz_0_zz[i] * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_0_xx, g_y_0_0_0_yz_yy_0_xy, g_y_0_0_0_yz_yy_0_xz, g_y_0_0_0_yz_yy_0_yy, g_y_0_0_0_yz_yy_0_yz, g_y_0_0_0_yz_yy_0_zz, g_yyz_yy_0_xx, g_yyz_yy_0_xy, g_yyz_yy_0_xz, g_yyz_yy_0_yy, g_yyz_yy_0_yz, g_yyz_yy_0_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_0_xx[i] = -g_z_yy_0_xx[i] + 2.0 * g_yyz_yy_0_xx[i] * a_exp;

        g_y_0_0_0_yz_yy_0_xy[i] = -g_z_yy_0_xy[i] + 2.0 * g_yyz_yy_0_xy[i] * a_exp;

        g_y_0_0_0_yz_yy_0_xz[i] = -g_z_yy_0_xz[i] + 2.0 * g_yyz_yy_0_xz[i] * a_exp;

        g_y_0_0_0_yz_yy_0_yy[i] = -g_z_yy_0_yy[i] + 2.0 * g_yyz_yy_0_yy[i] * a_exp;

        g_y_0_0_0_yz_yy_0_yz[i] = -g_z_yy_0_yz[i] + 2.0 * g_yyz_yy_0_yz[i] * a_exp;

        g_y_0_0_0_yz_yy_0_zz[i] = -g_z_yy_0_zz[i] + 2.0 * g_yyz_yy_0_zz[i] * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_0_xx, g_y_0_0_0_yz_yz_0_xy, g_y_0_0_0_yz_yz_0_xz, g_y_0_0_0_yz_yz_0_yy, g_y_0_0_0_yz_yz_0_yz, g_y_0_0_0_yz_yz_0_zz, g_yyz_yz_0_xx, g_yyz_yz_0_xy, g_yyz_yz_0_xz, g_yyz_yz_0_yy, g_yyz_yz_0_yz, g_yyz_yz_0_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_0_xx[i] = -g_z_yz_0_xx[i] + 2.0 * g_yyz_yz_0_xx[i] * a_exp;

        g_y_0_0_0_yz_yz_0_xy[i] = -g_z_yz_0_xy[i] + 2.0 * g_yyz_yz_0_xy[i] * a_exp;

        g_y_0_0_0_yz_yz_0_xz[i] = -g_z_yz_0_xz[i] + 2.0 * g_yyz_yz_0_xz[i] * a_exp;

        g_y_0_0_0_yz_yz_0_yy[i] = -g_z_yz_0_yy[i] + 2.0 * g_yyz_yz_0_yy[i] * a_exp;

        g_y_0_0_0_yz_yz_0_yz[i] = -g_z_yz_0_yz[i] + 2.0 * g_yyz_yz_0_yz[i] * a_exp;

        g_y_0_0_0_yz_yz_0_zz[i] = -g_z_yz_0_zz[i] + 2.0 * g_yyz_yz_0_zz[i] * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_0_xx, g_y_0_0_0_yz_zz_0_xy, g_y_0_0_0_yz_zz_0_xz, g_y_0_0_0_yz_zz_0_yy, g_y_0_0_0_yz_zz_0_yz, g_y_0_0_0_yz_zz_0_zz, g_yyz_zz_0_xx, g_yyz_zz_0_xy, g_yyz_zz_0_xz, g_yyz_zz_0_yy, g_yyz_zz_0_yz, g_yyz_zz_0_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_0_xx[i] = -g_z_zz_0_xx[i] + 2.0 * g_yyz_zz_0_xx[i] * a_exp;

        g_y_0_0_0_yz_zz_0_xy[i] = -g_z_zz_0_xy[i] + 2.0 * g_yyz_zz_0_xy[i] * a_exp;

        g_y_0_0_0_yz_zz_0_xz[i] = -g_z_zz_0_xz[i] + 2.0 * g_yyz_zz_0_xz[i] * a_exp;

        g_y_0_0_0_yz_zz_0_yy[i] = -g_z_zz_0_yy[i] + 2.0 * g_yyz_zz_0_yy[i] * a_exp;

        g_y_0_0_0_yz_zz_0_yz[i] = -g_z_zz_0_yz[i] + 2.0 * g_yyz_zz_0_yz[i] * a_exp;

        g_y_0_0_0_yz_zz_0_zz[i] = -g_z_zz_0_zz[i] + 2.0 * g_yyz_zz_0_zz[i] * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_0_xx, g_y_0_0_0_zz_xx_0_xy, g_y_0_0_0_zz_xx_0_xz, g_y_0_0_0_zz_xx_0_yy, g_y_0_0_0_zz_xx_0_yz, g_y_0_0_0_zz_xx_0_zz, g_yzz_xx_0_xx, g_yzz_xx_0_xy, g_yzz_xx_0_xz, g_yzz_xx_0_yy, g_yzz_xx_0_yz, g_yzz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_0_xx[i] = 2.0 * g_yzz_xx_0_xx[i] * a_exp;

        g_y_0_0_0_zz_xx_0_xy[i] = 2.0 * g_yzz_xx_0_xy[i] * a_exp;

        g_y_0_0_0_zz_xx_0_xz[i] = 2.0 * g_yzz_xx_0_xz[i] * a_exp;

        g_y_0_0_0_zz_xx_0_yy[i] = 2.0 * g_yzz_xx_0_yy[i] * a_exp;

        g_y_0_0_0_zz_xx_0_yz[i] = 2.0 * g_yzz_xx_0_yz[i] * a_exp;

        g_y_0_0_0_zz_xx_0_zz[i] = 2.0 * g_yzz_xx_0_zz[i] * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_0_xx, g_y_0_0_0_zz_xy_0_xy, g_y_0_0_0_zz_xy_0_xz, g_y_0_0_0_zz_xy_0_yy, g_y_0_0_0_zz_xy_0_yz, g_y_0_0_0_zz_xy_0_zz, g_yzz_xy_0_xx, g_yzz_xy_0_xy, g_yzz_xy_0_xz, g_yzz_xy_0_yy, g_yzz_xy_0_yz, g_yzz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_0_xx[i] = 2.0 * g_yzz_xy_0_xx[i] * a_exp;

        g_y_0_0_0_zz_xy_0_xy[i] = 2.0 * g_yzz_xy_0_xy[i] * a_exp;

        g_y_0_0_0_zz_xy_0_xz[i] = 2.0 * g_yzz_xy_0_xz[i] * a_exp;

        g_y_0_0_0_zz_xy_0_yy[i] = 2.0 * g_yzz_xy_0_yy[i] * a_exp;

        g_y_0_0_0_zz_xy_0_yz[i] = 2.0 * g_yzz_xy_0_yz[i] * a_exp;

        g_y_0_0_0_zz_xy_0_zz[i] = 2.0 * g_yzz_xy_0_zz[i] * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_0_xx, g_y_0_0_0_zz_xz_0_xy, g_y_0_0_0_zz_xz_0_xz, g_y_0_0_0_zz_xz_0_yy, g_y_0_0_0_zz_xz_0_yz, g_y_0_0_0_zz_xz_0_zz, g_yzz_xz_0_xx, g_yzz_xz_0_xy, g_yzz_xz_0_xz, g_yzz_xz_0_yy, g_yzz_xz_0_yz, g_yzz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_0_xx[i] = 2.0 * g_yzz_xz_0_xx[i] * a_exp;

        g_y_0_0_0_zz_xz_0_xy[i] = 2.0 * g_yzz_xz_0_xy[i] * a_exp;

        g_y_0_0_0_zz_xz_0_xz[i] = 2.0 * g_yzz_xz_0_xz[i] * a_exp;

        g_y_0_0_0_zz_xz_0_yy[i] = 2.0 * g_yzz_xz_0_yy[i] * a_exp;

        g_y_0_0_0_zz_xz_0_yz[i] = 2.0 * g_yzz_xz_0_yz[i] * a_exp;

        g_y_0_0_0_zz_xz_0_zz[i] = 2.0 * g_yzz_xz_0_zz[i] * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_0_xx, g_y_0_0_0_zz_yy_0_xy, g_y_0_0_0_zz_yy_0_xz, g_y_0_0_0_zz_yy_0_yy, g_y_0_0_0_zz_yy_0_yz, g_y_0_0_0_zz_yy_0_zz, g_yzz_yy_0_xx, g_yzz_yy_0_xy, g_yzz_yy_0_xz, g_yzz_yy_0_yy, g_yzz_yy_0_yz, g_yzz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_0_xx[i] = 2.0 * g_yzz_yy_0_xx[i] * a_exp;

        g_y_0_0_0_zz_yy_0_xy[i] = 2.0 * g_yzz_yy_0_xy[i] * a_exp;

        g_y_0_0_0_zz_yy_0_xz[i] = 2.0 * g_yzz_yy_0_xz[i] * a_exp;

        g_y_0_0_0_zz_yy_0_yy[i] = 2.0 * g_yzz_yy_0_yy[i] * a_exp;

        g_y_0_0_0_zz_yy_0_yz[i] = 2.0 * g_yzz_yy_0_yz[i] * a_exp;

        g_y_0_0_0_zz_yy_0_zz[i] = 2.0 * g_yzz_yy_0_zz[i] * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_0_xx, g_y_0_0_0_zz_yz_0_xy, g_y_0_0_0_zz_yz_0_xz, g_y_0_0_0_zz_yz_0_yy, g_y_0_0_0_zz_yz_0_yz, g_y_0_0_0_zz_yz_0_zz, g_yzz_yz_0_xx, g_yzz_yz_0_xy, g_yzz_yz_0_xz, g_yzz_yz_0_yy, g_yzz_yz_0_yz, g_yzz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_0_xx[i] = 2.0 * g_yzz_yz_0_xx[i] * a_exp;

        g_y_0_0_0_zz_yz_0_xy[i] = 2.0 * g_yzz_yz_0_xy[i] * a_exp;

        g_y_0_0_0_zz_yz_0_xz[i] = 2.0 * g_yzz_yz_0_xz[i] * a_exp;

        g_y_0_0_0_zz_yz_0_yy[i] = 2.0 * g_yzz_yz_0_yy[i] * a_exp;

        g_y_0_0_0_zz_yz_0_yz[i] = 2.0 * g_yzz_yz_0_yz[i] * a_exp;

        g_y_0_0_0_zz_yz_0_zz[i] = 2.0 * g_yzz_yz_0_zz[i] * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_0_xx, g_y_0_0_0_zz_zz_0_xy, g_y_0_0_0_zz_zz_0_xz, g_y_0_0_0_zz_zz_0_yy, g_y_0_0_0_zz_zz_0_yz, g_y_0_0_0_zz_zz_0_zz, g_yzz_zz_0_xx, g_yzz_zz_0_xy, g_yzz_zz_0_xz, g_yzz_zz_0_yy, g_yzz_zz_0_yz, g_yzz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_0_xx[i] = 2.0 * g_yzz_zz_0_xx[i] * a_exp;

        g_y_0_0_0_zz_zz_0_xy[i] = 2.0 * g_yzz_zz_0_xy[i] * a_exp;

        g_y_0_0_0_zz_zz_0_xz[i] = 2.0 * g_yzz_zz_0_xz[i] * a_exp;

        g_y_0_0_0_zz_zz_0_yy[i] = 2.0 * g_yzz_zz_0_yy[i] * a_exp;

        g_y_0_0_0_zz_zz_0_yz[i] = 2.0 * g_yzz_zz_0_yz[i] * a_exp;

        g_y_0_0_0_zz_zz_0_zz[i] = 2.0 * g_yzz_zz_0_zz[i] * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_xxz_xx_0_xx, g_xxz_xx_0_xy, g_xxz_xx_0_xz, g_xxz_xx_0_yy, g_xxz_xx_0_yz, g_xxz_xx_0_zz, g_z_0_0_0_xx_xx_0_xx, g_z_0_0_0_xx_xx_0_xy, g_z_0_0_0_xx_xx_0_xz, g_z_0_0_0_xx_xx_0_yy, g_z_0_0_0_xx_xx_0_yz, g_z_0_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_0_xx[i] = 2.0 * g_xxz_xx_0_xx[i] * a_exp;

        g_z_0_0_0_xx_xx_0_xy[i] = 2.0 * g_xxz_xx_0_xy[i] * a_exp;

        g_z_0_0_0_xx_xx_0_xz[i] = 2.0 * g_xxz_xx_0_xz[i] * a_exp;

        g_z_0_0_0_xx_xx_0_yy[i] = 2.0 * g_xxz_xx_0_yy[i] * a_exp;

        g_z_0_0_0_xx_xx_0_yz[i] = 2.0 * g_xxz_xx_0_yz[i] * a_exp;

        g_z_0_0_0_xx_xx_0_zz[i] = 2.0 * g_xxz_xx_0_zz[i] * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_xxz_xy_0_xx, g_xxz_xy_0_xy, g_xxz_xy_0_xz, g_xxz_xy_0_yy, g_xxz_xy_0_yz, g_xxz_xy_0_zz, g_z_0_0_0_xx_xy_0_xx, g_z_0_0_0_xx_xy_0_xy, g_z_0_0_0_xx_xy_0_xz, g_z_0_0_0_xx_xy_0_yy, g_z_0_0_0_xx_xy_0_yz, g_z_0_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_0_xx[i] = 2.0 * g_xxz_xy_0_xx[i] * a_exp;

        g_z_0_0_0_xx_xy_0_xy[i] = 2.0 * g_xxz_xy_0_xy[i] * a_exp;

        g_z_0_0_0_xx_xy_0_xz[i] = 2.0 * g_xxz_xy_0_xz[i] * a_exp;

        g_z_0_0_0_xx_xy_0_yy[i] = 2.0 * g_xxz_xy_0_yy[i] * a_exp;

        g_z_0_0_0_xx_xy_0_yz[i] = 2.0 * g_xxz_xy_0_yz[i] * a_exp;

        g_z_0_0_0_xx_xy_0_zz[i] = 2.0 * g_xxz_xy_0_zz[i] * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_xxz_xz_0_xx, g_xxz_xz_0_xy, g_xxz_xz_0_xz, g_xxz_xz_0_yy, g_xxz_xz_0_yz, g_xxz_xz_0_zz, g_z_0_0_0_xx_xz_0_xx, g_z_0_0_0_xx_xz_0_xy, g_z_0_0_0_xx_xz_0_xz, g_z_0_0_0_xx_xz_0_yy, g_z_0_0_0_xx_xz_0_yz, g_z_0_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_0_xx[i] = 2.0 * g_xxz_xz_0_xx[i] * a_exp;

        g_z_0_0_0_xx_xz_0_xy[i] = 2.0 * g_xxz_xz_0_xy[i] * a_exp;

        g_z_0_0_0_xx_xz_0_xz[i] = 2.0 * g_xxz_xz_0_xz[i] * a_exp;

        g_z_0_0_0_xx_xz_0_yy[i] = 2.0 * g_xxz_xz_0_yy[i] * a_exp;

        g_z_0_0_0_xx_xz_0_yz[i] = 2.0 * g_xxz_xz_0_yz[i] * a_exp;

        g_z_0_0_0_xx_xz_0_zz[i] = 2.0 * g_xxz_xz_0_zz[i] * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_xxz_yy_0_xx, g_xxz_yy_0_xy, g_xxz_yy_0_xz, g_xxz_yy_0_yy, g_xxz_yy_0_yz, g_xxz_yy_0_zz, g_z_0_0_0_xx_yy_0_xx, g_z_0_0_0_xx_yy_0_xy, g_z_0_0_0_xx_yy_0_xz, g_z_0_0_0_xx_yy_0_yy, g_z_0_0_0_xx_yy_0_yz, g_z_0_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_0_xx[i] = 2.0 * g_xxz_yy_0_xx[i] * a_exp;

        g_z_0_0_0_xx_yy_0_xy[i] = 2.0 * g_xxz_yy_0_xy[i] * a_exp;

        g_z_0_0_0_xx_yy_0_xz[i] = 2.0 * g_xxz_yy_0_xz[i] * a_exp;

        g_z_0_0_0_xx_yy_0_yy[i] = 2.0 * g_xxz_yy_0_yy[i] * a_exp;

        g_z_0_0_0_xx_yy_0_yz[i] = 2.0 * g_xxz_yy_0_yz[i] * a_exp;

        g_z_0_0_0_xx_yy_0_zz[i] = 2.0 * g_xxz_yy_0_zz[i] * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_xxz_yz_0_xx, g_xxz_yz_0_xy, g_xxz_yz_0_xz, g_xxz_yz_0_yy, g_xxz_yz_0_yz, g_xxz_yz_0_zz, g_z_0_0_0_xx_yz_0_xx, g_z_0_0_0_xx_yz_0_xy, g_z_0_0_0_xx_yz_0_xz, g_z_0_0_0_xx_yz_0_yy, g_z_0_0_0_xx_yz_0_yz, g_z_0_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_0_xx[i] = 2.0 * g_xxz_yz_0_xx[i] * a_exp;

        g_z_0_0_0_xx_yz_0_xy[i] = 2.0 * g_xxz_yz_0_xy[i] * a_exp;

        g_z_0_0_0_xx_yz_0_xz[i] = 2.0 * g_xxz_yz_0_xz[i] * a_exp;

        g_z_0_0_0_xx_yz_0_yy[i] = 2.0 * g_xxz_yz_0_yy[i] * a_exp;

        g_z_0_0_0_xx_yz_0_yz[i] = 2.0 * g_xxz_yz_0_yz[i] * a_exp;

        g_z_0_0_0_xx_yz_0_zz[i] = 2.0 * g_xxz_yz_0_zz[i] * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_xxz_zz_0_xx, g_xxz_zz_0_xy, g_xxz_zz_0_xz, g_xxz_zz_0_yy, g_xxz_zz_0_yz, g_xxz_zz_0_zz, g_z_0_0_0_xx_zz_0_xx, g_z_0_0_0_xx_zz_0_xy, g_z_0_0_0_xx_zz_0_xz, g_z_0_0_0_xx_zz_0_yy, g_z_0_0_0_xx_zz_0_yz, g_z_0_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_0_xx[i] = 2.0 * g_xxz_zz_0_xx[i] * a_exp;

        g_z_0_0_0_xx_zz_0_xy[i] = 2.0 * g_xxz_zz_0_xy[i] * a_exp;

        g_z_0_0_0_xx_zz_0_xz[i] = 2.0 * g_xxz_zz_0_xz[i] * a_exp;

        g_z_0_0_0_xx_zz_0_yy[i] = 2.0 * g_xxz_zz_0_yy[i] * a_exp;

        g_z_0_0_0_xx_zz_0_yz[i] = 2.0 * g_xxz_zz_0_yz[i] * a_exp;

        g_z_0_0_0_xx_zz_0_zz[i] = 2.0 * g_xxz_zz_0_zz[i] * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_xyz_xx_0_xx, g_xyz_xx_0_xy, g_xyz_xx_0_xz, g_xyz_xx_0_yy, g_xyz_xx_0_yz, g_xyz_xx_0_zz, g_z_0_0_0_xy_xx_0_xx, g_z_0_0_0_xy_xx_0_xy, g_z_0_0_0_xy_xx_0_xz, g_z_0_0_0_xy_xx_0_yy, g_z_0_0_0_xy_xx_0_yz, g_z_0_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_0_xx[i] = 2.0 * g_xyz_xx_0_xx[i] * a_exp;

        g_z_0_0_0_xy_xx_0_xy[i] = 2.0 * g_xyz_xx_0_xy[i] * a_exp;

        g_z_0_0_0_xy_xx_0_xz[i] = 2.0 * g_xyz_xx_0_xz[i] * a_exp;

        g_z_0_0_0_xy_xx_0_yy[i] = 2.0 * g_xyz_xx_0_yy[i] * a_exp;

        g_z_0_0_0_xy_xx_0_yz[i] = 2.0 * g_xyz_xx_0_yz[i] * a_exp;

        g_z_0_0_0_xy_xx_0_zz[i] = 2.0 * g_xyz_xx_0_zz[i] * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_xyz_xy_0_xx, g_xyz_xy_0_xy, g_xyz_xy_0_xz, g_xyz_xy_0_yy, g_xyz_xy_0_yz, g_xyz_xy_0_zz, g_z_0_0_0_xy_xy_0_xx, g_z_0_0_0_xy_xy_0_xy, g_z_0_0_0_xy_xy_0_xz, g_z_0_0_0_xy_xy_0_yy, g_z_0_0_0_xy_xy_0_yz, g_z_0_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_0_xx[i] = 2.0 * g_xyz_xy_0_xx[i] * a_exp;

        g_z_0_0_0_xy_xy_0_xy[i] = 2.0 * g_xyz_xy_0_xy[i] * a_exp;

        g_z_0_0_0_xy_xy_0_xz[i] = 2.0 * g_xyz_xy_0_xz[i] * a_exp;

        g_z_0_0_0_xy_xy_0_yy[i] = 2.0 * g_xyz_xy_0_yy[i] * a_exp;

        g_z_0_0_0_xy_xy_0_yz[i] = 2.0 * g_xyz_xy_0_yz[i] * a_exp;

        g_z_0_0_0_xy_xy_0_zz[i] = 2.0 * g_xyz_xy_0_zz[i] * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_xyz_xz_0_xx, g_xyz_xz_0_xy, g_xyz_xz_0_xz, g_xyz_xz_0_yy, g_xyz_xz_0_yz, g_xyz_xz_0_zz, g_z_0_0_0_xy_xz_0_xx, g_z_0_0_0_xy_xz_0_xy, g_z_0_0_0_xy_xz_0_xz, g_z_0_0_0_xy_xz_0_yy, g_z_0_0_0_xy_xz_0_yz, g_z_0_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_0_xx[i] = 2.0 * g_xyz_xz_0_xx[i] * a_exp;

        g_z_0_0_0_xy_xz_0_xy[i] = 2.0 * g_xyz_xz_0_xy[i] * a_exp;

        g_z_0_0_0_xy_xz_0_xz[i] = 2.0 * g_xyz_xz_0_xz[i] * a_exp;

        g_z_0_0_0_xy_xz_0_yy[i] = 2.0 * g_xyz_xz_0_yy[i] * a_exp;

        g_z_0_0_0_xy_xz_0_yz[i] = 2.0 * g_xyz_xz_0_yz[i] * a_exp;

        g_z_0_0_0_xy_xz_0_zz[i] = 2.0 * g_xyz_xz_0_zz[i] * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_xyz_yy_0_xx, g_xyz_yy_0_xy, g_xyz_yy_0_xz, g_xyz_yy_0_yy, g_xyz_yy_0_yz, g_xyz_yy_0_zz, g_z_0_0_0_xy_yy_0_xx, g_z_0_0_0_xy_yy_0_xy, g_z_0_0_0_xy_yy_0_xz, g_z_0_0_0_xy_yy_0_yy, g_z_0_0_0_xy_yy_0_yz, g_z_0_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_0_xx[i] = 2.0 * g_xyz_yy_0_xx[i] * a_exp;

        g_z_0_0_0_xy_yy_0_xy[i] = 2.0 * g_xyz_yy_0_xy[i] * a_exp;

        g_z_0_0_0_xy_yy_0_xz[i] = 2.0 * g_xyz_yy_0_xz[i] * a_exp;

        g_z_0_0_0_xy_yy_0_yy[i] = 2.0 * g_xyz_yy_0_yy[i] * a_exp;

        g_z_0_0_0_xy_yy_0_yz[i] = 2.0 * g_xyz_yy_0_yz[i] * a_exp;

        g_z_0_0_0_xy_yy_0_zz[i] = 2.0 * g_xyz_yy_0_zz[i] * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_xyz_yz_0_xx, g_xyz_yz_0_xy, g_xyz_yz_0_xz, g_xyz_yz_0_yy, g_xyz_yz_0_yz, g_xyz_yz_0_zz, g_z_0_0_0_xy_yz_0_xx, g_z_0_0_0_xy_yz_0_xy, g_z_0_0_0_xy_yz_0_xz, g_z_0_0_0_xy_yz_0_yy, g_z_0_0_0_xy_yz_0_yz, g_z_0_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_0_xx[i] = 2.0 * g_xyz_yz_0_xx[i] * a_exp;

        g_z_0_0_0_xy_yz_0_xy[i] = 2.0 * g_xyz_yz_0_xy[i] * a_exp;

        g_z_0_0_0_xy_yz_0_xz[i] = 2.0 * g_xyz_yz_0_xz[i] * a_exp;

        g_z_0_0_0_xy_yz_0_yy[i] = 2.0 * g_xyz_yz_0_yy[i] * a_exp;

        g_z_0_0_0_xy_yz_0_yz[i] = 2.0 * g_xyz_yz_0_yz[i] * a_exp;

        g_z_0_0_0_xy_yz_0_zz[i] = 2.0 * g_xyz_yz_0_zz[i] * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_xyz_zz_0_xx, g_xyz_zz_0_xy, g_xyz_zz_0_xz, g_xyz_zz_0_yy, g_xyz_zz_0_yz, g_xyz_zz_0_zz, g_z_0_0_0_xy_zz_0_xx, g_z_0_0_0_xy_zz_0_xy, g_z_0_0_0_xy_zz_0_xz, g_z_0_0_0_xy_zz_0_yy, g_z_0_0_0_xy_zz_0_yz, g_z_0_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_0_xx[i] = 2.0 * g_xyz_zz_0_xx[i] * a_exp;

        g_z_0_0_0_xy_zz_0_xy[i] = 2.0 * g_xyz_zz_0_xy[i] * a_exp;

        g_z_0_0_0_xy_zz_0_xz[i] = 2.0 * g_xyz_zz_0_xz[i] * a_exp;

        g_z_0_0_0_xy_zz_0_yy[i] = 2.0 * g_xyz_zz_0_yy[i] * a_exp;

        g_z_0_0_0_xy_zz_0_yz[i] = 2.0 * g_xyz_zz_0_yz[i] * a_exp;

        g_z_0_0_0_xy_zz_0_zz[i] = 2.0 * g_xyz_zz_0_zz[i] * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_xzz_xx_0_xx, g_xzz_xx_0_xy, g_xzz_xx_0_xz, g_xzz_xx_0_yy, g_xzz_xx_0_yz, g_xzz_xx_0_zz, g_z_0_0_0_xz_xx_0_xx, g_z_0_0_0_xz_xx_0_xy, g_z_0_0_0_xz_xx_0_xz, g_z_0_0_0_xz_xx_0_yy, g_z_0_0_0_xz_xx_0_yz, g_z_0_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_0_xx[i] = -g_x_xx_0_xx[i] + 2.0 * g_xzz_xx_0_xx[i] * a_exp;

        g_z_0_0_0_xz_xx_0_xy[i] = -g_x_xx_0_xy[i] + 2.0 * g_xzz_xx_0_xy[i] * a_exp;

        g_z_0_0_0_xz_xx_0_xz[i] = -g_x_xx_0_xz[i] + 2.0 * g_xzz_xx_0_xz[i] * a_exp;

        g_z_0_0_0_xz_xx_0_yy[i] = -g_x_xx_0_yy[i] + 2.0 * g_xzz_xx_0_yy[i] * a_exp;

        g_z_0_0_0_xz_xx_0_yz[i] = -g_x_xx_0_yz[i] + 2.0 * g_xzz_xx_0_yz[i] * a_exp;

        g_z_0_0_0_xz_xx_0_zz[i] = -g_x_xx_0_zz[i] + 2.0 * g_xzz_xx_0_zz[i] * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_xzz_xy_0_xx, g_xzz_xy_0_xy, g_xzz_xy_0_xz, g_xzz_xy_0_yy, g_xzz_xy_0_yz, g_xzz_xy_0_zz, g_z_0_0_0_xz_xy_0_xx, g_z_0_0_0_xz_xy_0_xy, g_z_0_0_0_xz_xy_0_xz, g_z_0_0_0_xz_xy_0_yy, g_z_0_0_0_xz_xy_0_yz, g_z_0_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_0_xx[i] = -g_x_xy_0_xx[i] + 2.0 * g_xzz_xy_0_xx[i] * a_exp;

        g_z_0_0_0_xz_xy_0_xy[i] = -g_x_xy_0_xy[i] + 2.0 * g_xzz_xy_0_xy[i] * a_exp;

        g_z_0_0_0_xz_xy_0_xz[i] = -g_x_xy_0_xz[i] + 2.0 * g_xzz_xy_0_xz[i] * a_exp;

        g_z_0_0_0_xz_xy_0_yy[i] = -g_x_xy_0_yy[i] + 2.0 * g_xzz_xy_0_yy[i] * a_exp;

        g_z_0_0_0_xz_xy_0_yz[i] = -g_x_xy_0_yz[i] + 2.0 * g_xzz_xy_0_yz[i] * a_exp;

        g_z_0_0_0_xz_xy_0_zz[i] = -g_x_xy_0_zz[i] + 2.0 * g_xzz_xy_0_zz[i] * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_xzz_xz_0_xx, g_xzz_xz_0_xy, g_xzz_xz_0_xz, g_xzz_xz_0_yy, g_xzz_xz_0_yz, g_xzz_xz_0_zz, g_z_0_0_0_xz_xz_0_xx, g_z_0_0_0_xz_xz_0_xy, g_z_0_0_0_xz_xz_0_xz, g_z_0_0_0_xz_xz_0_yy, g_z_0_0_0_xz_xz_0_yz, g_z_0_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_0_xx[i] = -g_x_xz_0_xx[i] + 2.0 * g_xzz_xz_0_xx[i] * a_exp;

        g_z_0_0_0_xz_xz_0_xy[i] = -g_x_xz_0_xy[i] + 2.0 * g_xzz_xz_0_xy[i] * a_exp;

        g_z_0_0_0_xz_xz_0_xz[i] = -g_x_xz_0_xz[i] + 2.0 * g_xzz_xz_0_xz[i] * a_exp;

        g_z_0_0_0_xz_xz_0_yy[i] = -g_x_xz_0_yy[i] + 2.0 * g_xzz_xz_0_yy[i] * a_exp;

        g_z_0_0_0_xz_xz_0_yz[i] = -g_x_xz_0_yz[i] + 2.0 * g_xzz_xz_0_yz[i] * a_exp;

        g_z_0_0_0_xz_xz_0_zz[i] = -g_x_xz_0_zz[i] + 2.0 * g_xzz_xz_0_zz[i] * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_xzz_yy_0_xx, g_xzz_yy_0_xy, g_xzz_yy_0_xz, g_xzz_yy_0_yy, g_xzz_yy_0_yz, g_xzz_yy_0_zz, g_z_0_0_0_xz_yy_0_xx, g_z_0_0_0_xz_yy_0_xy, g_z_0_0_0_xz_yy_0_xz, g_z_0_0_0_xz_yy_0_yy, g_z_0_0_0_xz_yy_0_yz, g_z_0_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_0_xx[i] = -g_x_yy_0_xx[i] + 2.0 * g_xzz_yy_0_xx[i] * a_exp;

        g_z_0_0_0_xz_yy_0_xy[i] = -g_x_yy_0_xy[i] + 2.0 * g_xzz_yy_0_xy[i] * a_exp;

        g_z_0_0_0_xz_yy_0_xz[i] = -g_x_yy_0_xz[i] + 2.0 * g_xzz_yy_0_xz[i] * a_exp;

        g_z_0_0_0_xz_yy_0_yy[i] = -g_x_yy_0_yy[i] + 2.0 * g_xzz_yy_0_yy[i] * a_exp;

        g_z_0_0_0_xz_yy_0_yz[i] = -g_x_yy_0_yz[i] + 2.0 * g_xzz_yy_0_yz[i] * a_exp;

        g_z_0_0_0_xz_yy_0_zz[i] = -g_x_yy_0_zz[i] + 2.0 * g_xzz_yy_0_zz[i] * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_xzz_yz_0_xx, g_xzz_yz_0_xy, g_xzz_yz_0_xz, g_xzz_yz_0_yy, g_xzz_yz_0_yz, g_xzz_yz_0_zz, g_z_0_0_0_xz_yz_0_xx, g_z_0_0_0_xz_yz_0_xy, g_z_0_0_0_xz_yz_0_xz, g_z_0_0_0_xz_yz_0_yy, g_z_0_0_0_xz_yz_0_yz, g_z_0_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_0_xx[i] = -g_x_yz_0_xx[i] + 2.0 * g_xzz_yz_0_xx[i] * a_exp;

        g_z_0_0_0_xz_yz_0_xy[i] = -g_x_yz_0_xy[i] + 2.0 * g_xzz_yz_0_xy[i] * a_exp;

        g_z_0_0_0_xz_yz_0_xz[i] = -g_x_yz_0_xz[i] + 2.0 * g_xzz_yz_0_xz[i] * a_exp;

        g_z_0_0_0_xz_yz_0_yy[i] = -g_x_yz_0_yy[i] + 2.0 * g_xzz_yz_0_yy[i] * a_exp;

        g_z_0_0_0_xz_yz_0_yz[i] = -g_x_yz_0_yz[i] + 2.0 * g_xzz_yz_0_yz[i] * a_exp;

        g_z_0_0_0_xz_yz_0_zz[i] = -g_x_yz_0_zz[i] + 2.0 * g_xzz_yz_0_zz[i] * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_xzz_zz_0_xx, g_xzz_zz_0_xy, g_xzz_zz_0_xz, g_xzz_zz_0_yy, g_xzz_zz_0_yz, g_xzz_zz_0_zz, g_z_0_0_0_xz_zz_0_xx, g_z_0_0_0_xz_zz_0_xy, g_z_0_0_0_xz_zz_0_xz, g_z_0_0_0_xz_zz_0_yy, g_z_0_0_0_xz_zz_0_yz, g_z_0_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_0_xx[i] = -g_x_zz_0_xx[i] + 2.0 * g_xzz_zz_0_xx[i] * a_exp;

        g_z_0_0_0_xz_zz_0_xy[i] = -g_x_zz_0_xy[i] + 2.0 * g_xzz_zz_0_xy[i] * a_exp;

        g_z_0_0_0_xz_zz_0_xz[i] = -g_x_zz_0_xz[i] + 2.0 * g_xzz_zz_0_xz[i] * a_exp;

        g_z_0_0_0_xz_zz_0_yy[i] = -g_x_zz_0_yy[i] + 2.0 * g_xzz_zz_0_yy[i] * a_exp;

        g_z_0_0_0_xz_zz_0_yz[i] = -g_x_zz_0_yz[i] + 2.0 * g_xzz_zz_0_yz[i] * a_exp;

        g_z_0_0_0_xz_zz_0_zz[i] = -g_x_zz_0_zz[i] + 2.0 * g_xzz_zz_0_zz[i] * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_yyz_xx_0_xx, g_yyz_xx_0_xy, g_yyz_xx_0_xz, g_yyz_xx_0_yy, g_yyz_xx_0_yz, g_yyz_xx_0_zz, g_z_0_0_0_yy_xx_0_xx, g_z_0_0_0_yy_xx_0_xy, g_z_0_0_0_yy_xx_0_xz, g_z_0_0_0_yy_xx_0_yy, g_z_0_0_0_yy_xx_0_yz, g_z_0_0_0_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_0_xx[i] = 2.0 * g_yyz_xx_0_xx[i] * a_exp;

        g_z_0_0_0_yy_xx_0_xy[i] = 2.0 * g_yyz_xx_0_xy[i] * a_exp;

        g_z_0_0_0_yy_xx_0_xz[i] = 2.0 * g_yyz_xx_0_xz[i] * a_exp;

        g_z_0_0_0_yy_xx_0_yy[i] = 2.0 * g_yyz_xx_0_yy[i] * a_exp;

        g_z_0_0_0_yy_xx_0_yz[i] = 2.0 * g_yyz_xx_0_yz[i] * a_exp;

        g_z_0_0_0_yy_xx_0_zz[i] = 2.0 * g_yyz_xx_0_zz[i] * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_yyz_xy_0_xx, g_yyz_xy_0_xy, g_yyz_xy_0_xz, g_yyz_xy_0_yy, g_yyz_xy_0_yz, g_yyz_xy_0_zz, g_z_0_0_0_yy_xy_0_xx, g_z_0_0_0_yy_xy_0_xy, g_z_0_0_0_yy_xy_0_xz, g_z_0_0_0_yy_xy_0_yy, g_z_0_0_0_yy_xy_0_yz, g_z_0_0_0_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_0_xx[i] = 2.0 * g_yyz_xy_0_xx[i] * a_exp;

        g_z_0_0_0_yy_xy_0_xy[i] = 2.0 * g_yyz_xy_0_xy[i] * a_exp;

        g_z_0_0_0_yy_xy_0_xz[i] = 2.0 * g_yyz_xy_0_xz[i] * a_exp;

        g_z_0_0_0_yy_xy_0_yy[i] = 2.0 * g_yyz_xy_0_yy[i] * a_exp;

        g_z_0_0_0_yy_xy_0_yz[i] = 2.0 * g_yyz_xy_0_yz[i] * a_exp;

        g_z_0_0_0_yy_xy_0_zz[i] = 2.0 * g_yyz_xy_0_zz[i] * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_yyz_xz_0_xx, g_yyz_xz_0_xy, g_yyz_xz_0_xz, g_yyz_xz_0_yy, g_yyz_xz_0_yz, g_yyz_xz_0_zz, g_z_0_0_0_yy_xz_0_xx, g_z_0_0_0_yy_xz_0_xy, g_z_0_0_0_yy_xz_0_xz, g_z_0_0_0_yy_xz_0_yy, g_z_0_0_0_yy_xz_0_yz, g_z_0_0_0_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_0_xx[i] = 2.0 * g_yyz_xz_0_xx[i] * a_exp;

        g_z_0_0_0_yy_xz_0_xy[i] = 2.0 * g_yyz_xz_0_xy[i] * a_exp;

        g_z_0_0_0_yy_xz_0_xz[i] = 2.0 * g_yyz_xz_0_xz[i] * a_exp;

        g_z_0_0_0_yy_xz_0_yy[i] = 2.0 * g_yyz_xz_0_yy[i] * a_exp;

        g_z_0_0_0_yy_xz_0_yz[i] = 2.0 * g_yyz_xz_0_yz[i] * a_exp;

        g_z_0_0_0_yy_xz_0_zz[i] = 2.0 * g_yyz_xz_0_zz[i] * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_yyz_yy_0_xx, g_yyz_yy_0_xy, g_yyz_yy_0_xz, g_yyz_yy_0_yy, g_yyz_yy_0_yz, g_yyz_yy_0_zz, g_z_0_0_0_yy_yy_0_xx, g_z_0_0_0_yy_yy_0_xy, g_z_0_0_0_yy_yy_0_xz, g_z_0_0_0_yy_yy_0_yy, g_z_0_0_0_yy_yy_0_yz, g_z_0_0_0_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_0_xx[i] = 2.0 * g_yyz_yy_0_xx[i] * a_exp;

        g_z_0_0_0_yy_yy_0_xy[i] = 2.0 * g_yyz_yy_0_xy[i] * a_exp;

        g_z_0_0_0_yy_yy_0_xz[i] = 2.0 * g_yyz_yy_0_xz[i] * a_exp;

        g_z_0_0_0_yy_yy_0_yy[i] = 2.0 * g_yyz_yy_0_yy[i] * a_exp;

        g_z_0_0_0_yy_yy_0_yz[i] = 2.0 * g_yyz_yy_0_yz[i] * a_exp;

        g_z_0_0_0_yy_yy_0_zz[i] = 2.0 * g_yyz_yy_0_zz[i] * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_yyz_yz_0_xx, g_yyz_yz_0_xy, g_yyz_yz_0_xz, g_yyz_yz_0_yy, g_yyz_yz_0_yz, g_yyz_yz_0_zz, g_z_0_0_0_yy_yz_0_xx, g_z_0_0_0_yy_yz_0_xy, g_z_0_0_0_yy_yz_0_xz, g_z_0_0_0_yy_yz_0_yy, g_z_0_0_0_yy_yz_0_yz, g_z_0_0_0_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_0_xx[i] = 2.0 * g_yyz_yz_0_xx[i] * a_exp;

        g_z_0_0_0_yy_yz_0_xy[i] = 2.0 * g_yyz_yz_0_xy[i] * a_exp;

        g_z_0_0_0_yy_yz_0_xz[i] = 2.0 * g_yyz_yz_0_xz[i] * a_exp;

        g_z_0_0_0_yy_yz_0_yy[i] = 2.0 * g_yyz_yz_0_yy[i] * a_exp;

        g_z_0_0_0_yy_yz_0_yz[i] = 2.0 * g_yyz_yz_0_yz[i] * a_exp;

        g_z_0_0_0_yy_yz_0_zz[i] = 2.0 * g_yyz_yz_0_zz[i] * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_yyz_zz_0_xx, g_yyz_zz_0_xy, g_yyz_zz_0_xz, g_yyz_zz_0_yy, g_yyz_zz_0_yz, g_yyz_zz_0_zz, g_z_0_0_0_yy_zz_0_xx, g_z_0_0_0_yy_zz_0_xy, g_z_0_0_0_yy_zz_0_xz, g_z_0_0_0_yy_zz_0_yy, g_z_0_0_0_yy_zz_0_yz, g_z_0_0_0_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_0_xx[i] = 2.0 * g_yyz_zz_0_xx[i] * a_exp;

        g_z_0_0_0_yy_zz_0_xy[i] = 2.0 * g_yyz_zz_0_xy[i] * a_exp;

        g_z_0_0_0_yy_zz_0_xz[i] = 2.0 * g_yyz_zz_0_xz[i] * a_exp;

        g_z_0_0_0_yy_zz_0_yy[i] = 2.0 * g_yyz_zz_0_yy[i] * a_exp;

        g_z_0_0_0_yy_zz_0_yz[i] = 2.0 * g_yyz_zz_0_yz[i] * a_exp;

        g_z_0_0_0_yy_zz_0_zz[i] = 2.0 * g_yyz_zz_0_zz[i] * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz, g_yzz_xx_0_xx, g_yzz_xx_0_xy, g_yzz_xx_0_xz, g_yzz_xx_0_yy, g_yzz_xx_0_yz, g_yzz_xx_0_zz, g_z_0_0_0_yz_xx_0_xx, g_z_0_0_0_yz_xx_0_xy, g_z_0_0_0_yz_xx_0_xz, g_z_0_0_0_yz_xx_0_yy, g_z_0_0_0_yz_xx_0_yz, g_z_0_0_0_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_0_xx[i] = -g_y_xx_0_xx[i] + 2.0 * g_yzz_xx_0_xx[i] * a_exp;

        g_z_0_0_0_yz_xx_0_xy[i] = -g_y_xx_0_xy[i] + 2.0 * g_yzz_xx_0_xy[i] * a_exp;

        g_z_0_0_0_yz_xx_0_xz[i] = -g_y_xx_0_xz[i] + 2.0 * g_yzz_xx_0_xz[i] * a_exp;

        g_z_0_0_0_yz_xx_0_yy[i] = -g_y_xx_0_yy[i] + 2.0 * g_yzz_xx_0_yy[i] * a_exp;

        g_z_0_0_0_yz_xx_0_yz[i] = -g_y_xx_0_yz[i] + 2.0 * g_yzz_xx_0_yz[i] * a_exp;

        g_z_0_0_0_yz_xx_0_zz[i] = -g_y_xx_0_zz[i] + 2.0 * g_yzz_xx_0_zz[i] * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz, g_yzz_xy_0_xx, g_yzz_xy_0_xy, g_yzz_xy_0_xz, g_yzz_xy_0_yy, g_yzz_xy_0_yz, g_yzz_xy_0_zz, g_z_0_0_0_yz_xy_0_xx, g_z_0_0_0_yz_xy_0_xy, g_z_0_0_0_yz_xy_0_xz, g_z_0_0_0_yz_xy_0_yy, g_z_0_0_0_yz_xy_0_yz, g_z_0_0_0_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_0_xx[i] = -g_y_xy_0_xx[i] + 2.0 * g_yzz_xy_0_xx[i] * a_exp;

        g_z_0_0_0_yz_xy_0_xy[i] = -g_y_xy_0_xy[i] + 2.0 * g_yzz_xy_0_xy[i] * a_exp;

        g_z_0_0_0_yz_xy_0_xz[i] = -g_y_xy_0_xz[i] + 2.0 * g_yzz_xy_0_xz[i] * a_exp;

        g_z_0_0_0_yz_xy_0_yy[i] = -g_y_xy_0_yy[i] + 2.0 * g_yzz_xy_0_yy[i] * a_exp;

        g_z_0_0_0_yz_xy_0_yz[i] = -g_y_xy_0_yz[i] + 2.0 * g_yzz_xy_0_yz[i] * a_exp;

        g_z_0_0_0_yz_xy_0_zz[i] = -g_y_xy_0_zz[i] + 2.0 * g_yzz_xy_0_zz[i] * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz, g_yzz_xz_0_xx, g_yzz_xz_0_xy, g_yzz_xz_0_xz, g_yzz_xz_0_yy, g_yzz_xz_0_yz, g_yzz_xz_0_zz, g_z_0_0_0_yz_xz_0_xx, g_z_0_0_0_yz_xz_0_xy, g_z_0_0_0_yz_xz_0_xz, g_z_0_0_0_yz_xz_0_yy, g_z_0_0_0_yz_xz_0_yz, g_z_0_0_0_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_0_xx[i] = -g_y_xz_0_xx[i] + 2.0 * g_yzz_xz_0_xx[i] * a_exp;

        g_z_0_0_0_yz_xz_0_xy[i] = -g_y_xz_0_xy[i] + 2.0 * g_yzz_xz_0_xy[i] * a_exp;

        g_z_0_0_0_yz_xz_0_xz[i] = -g_y_xz_0_xz[i] + 2.0 * g_yzz_xz_0_xz[i] * a_exp;

        g_z_0_0_0_yz_xz_0_yy[i] = -g_y_xz_0_yy[i] + 2.0 * g_yzz_xz_0_yy[i] * a_exp;

        g_z_0_0_0_yz_xz_0_yz[i] = -g_y_xz_0_yz[i] + 2.0 * g_yzz_xz_0_yz[i] * a_exp;

        g_z_0_0_0_yz_xz_0_zz[i] = -g_y_xz_0_zz[i] + 2.0 * g_yzz_xz_0_zz[i] * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz, g_yzz_yy_0_xx, g_yzz_yy_0_xy, g_yzz_yy_0_xz, g_yzz_yy_0_yy, g_yzz_yy_0_yz, g_yzz_yy_0_zz, g_z_0_0_0_yz_yy_0_xx, g_z_0_0_0_yz_yy_0_xy, g_z_0_0_0_yz_yy_0_xz, g_z_0_0_0_yz_yy_0_yy, g_z_0_0_0_yz_yy_0_yz, g_z_0_0_0_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_0_xx[i] = -g_y_yy_0_xx[i] + 2.0 * g_yzz_yy_0_xx[i] * a_exp;

        g_z_0_0_0_yz_yy_0_xy[i] = -g_y_yy_0_xy[i] + 2.0 * g_yzz_yy_0_xy[i] * a_exp;

        g_z_0_0_0_yz_yy_0_xz[i] = -g_y_yy_0_xz[i] + 2.0 * g_yzz_yy_0_xz[i] * a_exp;

        g_z_0_0_0_yz_yy_0_yy[i] = -g_y_yy_0_yy[i] + 2.0 * g_yzz_yy_0_yy[i] * a_exp;

        g_z_0_0_0_yz_yy_0_yz[i] = -g_y_yy_0_yz[i] + 2.0 * g_yzz_yy_0_yz[i] * a_exp;

        g_z_0_0_0_yz_yy_0_zz[i] = -g_y_yy_0_zz[i] + 2.0 * g_yzz_yy_0_zz[i] * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz, g_yzz_yz_0_xx, g_yzz_yz_0_xy, g_yzz_yz_0_xz, g_yzz_yz_0_yy, g_yzz_yz_0_yz, g_yzz_yz_0_zz, g_z_0_0_0_yz_yz_0_xx, g_z_0_0_0_yz_yz_0_xy, g_z_0_0_0_yz_yz_0_xz, g_z_0_0_0_yz_yz_0_yy, g_z_0_0_0_yz_yz_0_yz, g_z_0_0_0_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_0_xx[i] = -g_y_yz_0_xx[i] + 2.0 * g_yzz_yz_0_xx[i] * a_exp;

        g_z_0_0_0_yz_yz_0_xy[i] = -g_y_yz_0_xy[i] + 2.0 * g_yzz_yz_0_xy[i] * a_exp;

        g_z_0_0_0_yz_yz_0_xz[i] = -g_y_yz_0_xz[i] + 2.0 * g_yzz_yz_0_xz[i] * a_exp;

        g_z_0_0_0_yz_yz_0_yy[i] = -g_y_yz_0_yy[i] + 2.0 * g_yzz_yz_0_yy[i] * a_exp;

        g_z_0_0_0_yz_yz_0_yz[i] = -g_y_yz_0_yz[i] + 2.0 * g_yzz_yz_0_yz[i] * a_exp;

        g_z_0_0_0_yz_yz_0_zz[i] = -g_y_yz_0_zz[i] + 2.0 * g_yzz_yz_0_zz[i] * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz, g_yzz_zz_0_xx, g_yzz_zz_0_xy, g_yzz_zz_0_xz, g_yzz_zz_0_yy, g_yzz_zz_0_yz, g_yzz_zz_0_zz, g_z_0_0_0_yz_zz_0_xx, g_z_0_0_0_yz_zz_0_xy, g_z_0_0_0_yz_zz_0_xz, g_z_0_0_0_yz_zz_0_yy, g_z_0_0_0_yz_zz_0_yz, g_z_0_0_0_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_0_xx[i] = -g_y_zz_0_xx[i] + 2.0 * g_yzz_zz_0_xx[i] * a_exp;

        g_z_0_0_0_yz_zz_0_xy[i] = -g_y_zz_0_xy[i] + 2.0 * g_yzz_zz_0_xy[i] * a_exp;

        g_z_0_0_0_yz_zz_0_xz[i] = -g_y_zz_0_xz[i] + 2.0 * g_yzz_zz_0_xz[i] * a_exp;

        g_z_0_0_0_yz_zz_0_yy[i] = -g_y_zz_0_yy[i] + 2.0 * g_yzz_zz_0_yy[i] * a_exp;

        g_z_0_0_0_yz_zz_0_yz[i] = -g_y_zz_0_yz[i] + 2.0 * g_yzz_zz_0_yz[i] * a_exp;

        g_z_0_0_0_yz_zz_0_zz[i] = -g_y_zz_0_zz[i] + 2.0 * g_yzz_zz_0_zz[i] * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_0_xx, g_z_0_0_0_zz_xx_0_xy, g_z_0_0_0_zz_xx_0_xz, g_z_0_0_0_zz_xx_0_yy, g_z_0_0_0_zz_xx_0_yz, g_z_0_0_0_zz_xx_0_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz, g_zzz_xx_0_xx, g_zzz_xx_0_xy, g_zzz_xx_0_xz, g_zzz_xx_0_yy, g_zzz_xx_0_yz, g_zzz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_0_xx[i] = -2.0 * g_z_xx_0_xx[i] + 2.0 * g_zzz_xx_0_xx[i] * a_exp;

        g_z_0_0_0_zz_xx_0_xy[i] = -2.0 * g_z_xx_0_xy[i] + 2.0 * g_zzz_xx_0_xy[i] * a_exp;

        g_z_0_0_0_zz_xx_0_xz[i] = -2.0 * g_z_xx_0_xz[i] + 2.0 * g_zzz_xx_0_xz[i] * a_exp;

        g_z_0_0_0_zz_xx_0_yy[i] = -2.0 * g_z_xx_0_yy[i] + 2.0 * g_zzz_xx_0_yy[i] * a_exp;

        g_z_0_0_0_zz_xx_0_yz[i] = -2.0 * g_z_xx_0_yz[i] + 2.0 * g_zzz_xx_0_yz[i] * a_exp;

        g_z_0_0_0_zz_xx_0_zz[i] = -2.0 * g_z_xx_0_zz[i] + 2.0 * g_zzz_xx_0_zz[i] * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_0_xx, g_z_0_0_0_zz_xy_0_xy, g_z_0_0_0_zz_xy_0_xz, g_z_0_0_0_zz_xy_0_yy, g_z_0_0_0_zz_xy_0_yz, g_z_0_0_0_zz_xy_0_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz, g_zzz_xy_0_xx, g_zzz_xy_0_xy, g_zzz_xy_0_xz, g_zzz_xy_0_yy, g_zzz_xy_0_yz, g_zzz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_0_xx[i] = -2.0 * g_z_xy_0_xx[i] + 2.0 * g_zzz_xy_0_xx[i] * a_exp;

        g_z_0_0_0_zz_xy_0_xy[i] = -2.0 * g_z_xy_0_xy[i] + 2.0 * g_zzz_xy_0_xy[i] * a_exp;

        g_z_0_0_0_zz_xy_0_xz[i] = -2.0 * g_z_xy_0_xz[i] + 2.0 * g_zzz_xy_0_xz[i] * a_exp;

        g_z_0_0_0_zz_xy_0_yy[i] = -2.0 * g_z_xy_0_yy[i] + 2.0 * g_zzz_xy_0_yy[i] * a_exp;

        g_z_0_0_0_zz_xy_0_yz[i] = -2.0 * g_z_xy_0_yz[i] + 2.0 * g_zzz_xy_0_yz[i] * a_exp;

        g_z_0_0_0_zz_xy_0_zz[i] = -2.0 * g_z_xy_0_zz[i] + 2.0 * g_zzz_xy_0_zz[i] * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_0_xx, g_z_0_0_0_zz_xz_0_xy, g_z_0_0_0_zz_xz_0_xz, g_z_0_0_0_zz_xz_0_yy, g_z_0_0_0_zz_xz_0_yz, g_z_0_0_0_zz_xz_0_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz, g_zzz_xz_0_xx, g_zzz_xz_0_xy, g_zzz_xz_0_xz, g_zzz_xz_0_yy, g_zzz_xz_0_yz, g_zzz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_0_xx[i] = -2.0 * g_z_xz_0_xx[i] + 2.0 * g_zzz_xz_0_xx[i] * a_exp;

        g_z_0_0_0_zz_xz_0_xy[i] = -2.0 * g_z_xz_0_xy[i] + 2.0 * g_zzz_xz_0_xy[i] * a_exp;

        g_z_0_0_0_zz_xz_0_xz[i] = -2.0 * g_z_xz_0_xz[i] + 2.0 * g_zzz_xz_0_xz[i] * a_exp;

        g_z_0_0_0_zz_xz_0_yy[i] = -2.0 * g_z_xz_0_yy[i] + 2.0 * g_zzz_xz_0_yy[i] * a_exp;

        g_z_0_0_0_zz_xz_0_yz[i] = -2.0 * g_z_xz_0_yz[i] + 2.0 * g_zzz_xz_0_yz[i] * a_exp;

        g_z_0_0_0_zz_xz_0_zz[i] = -2.0 * g_z_xz_0_zz[i] + 2.0 * g_zzz_xz_0_zz[i] * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_0_xx, g_z_0_0_0_zz_yy_0_xy, g_z_0_0_0_zz_yy_0_xz, g_z_0_0_0_zz_yy_0_yy, g_z_0_0_0_zz_yy_0_yz, g_z_0_0_0_zz_yy_0_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz, g_zzz_yy_0_xx, g_zzz_yy_0_xy, g_zzz_yy_0_xz, g_zzz_yy_0_yy, g_zzz_yy_0_yz, g_zzz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_0_xx[i] = -2.0 * g_z_yy_0_xx[i] + 2.0 * g_zzz_yy_0_xx[i] * a_exp;

        g_z_0_0_0_zz_yy_0_xy[i] = -2.0 * g_z_yy_0_xy[i] + 2.0 * g_zzz_yy_0_xy[i] * a_exp;

        g_z_0_0_0_zz_yy_0_xz[i] = -2.0 * g_z_yy_0_xz[i] + 2.0 * g_zzz_yy_0_xz[i] * a_exp;

        g_z_0_0_0_zz_yy_0_yy[i] = -2.0 * g_z_yy_0_yy[i] + 2.0 * g_zzz_yy_0_yy[i] * a_exp;

        g_z_0_0_0_zz_yy_0_yz[i] = -2.0 * g_z_yy_0_yz[i] + 2.0 * g_zzz_yy_0_yz[i] * a_exp;

        g_z_0_0_0_zz_yy_0_zz[i] = -2.0 * g_z_yy_0_zz[i] + 2.0 * g_zzz_yy_0_zz[i] * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_0_xx, g_z_0_0_0_zz_yz_0_xy, g_z_0_0_0_zz_yz_0_xz, g_z_0_0_0_zz_yz_0_yy, g_z_0_0_0_zz_yz_0_yz, g_z_0_0_0_zz_yz_0_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz, g_zzz_yz_0_xx, g_zzz_yz_0_xy, g_zzz_yz_0_xz, g_zzz_yz_0_yy, g_zzz_yz_0_yz, g_zzz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_0_xx[i] = -2.0 * g_z_yz_0_xx[i] + 2.0 * g_zzz_yz_0_xx[i] * a_exp;

        g_z_0_0_0_zz_yz_0_xy[i] = -2.0 * g_z_yz_0_xy[i] + 2.0 * g_zzz_yz_0_xy[i] * a_exp;

        g_z_0_0_0_zz_yz_0_xz[i] = -2.0 * g_z_yz_0_xz[i] + 2.0 * g_zzz_yz_0_xz[i] * a_exp;

        g_z_0_0_0_zz_yz_0_yy[i] = -2.0 * g_z_yz_0_yy[i] + 2.0 * g_zzz_yz_0_yy[i] * a_exp;

        g_z_0_0_0_zz_yz_0_yz[i] = -2.0 * g_z_yz_0_yz[i] + 2.0 * g_zzz_yz_0_yz[i] * a_exp;

        g_z_0_0_0_zz_yz_0_zz[i] = -2.0 * g_z_yz_0_zz[i] + 2.0 * g_zzz_yz_0_zz[i] * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_0_xx, g_z_0_0_0_zz_zz_0_xy, g_z_0_0_0_zz_zz_0_xz, g_z_0_0_0_zz_zz_0_yy, g_z_0_0_0_zz_zz_0_yz, g_z_0_0_0_zz_zz_0_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz, g_zzz_zz_0_xx, g_zzz_zz_0_xy, g_zzz_zz_0_xz, g_zzz_zz_0_yy, g_zzz_zz_0_yz, g_zzz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_0_xx[i] = -2.0 * g_z_zz_0_xx[i] + 2.0 * g_zzz_zz_0_xx[i] * a_exp;

        g_z_0_0_0_zz_zz_0_xy[i] = -2.0 * g_z_zz_0_xy[i] + 2.0 * g_zzz_zz_0_xy[i] * a_exp;

        g_z_0_0_0_zz_zz_0_xz[i] = -2.0 * g_z_zz_0_xz[i] + 2.0 * g_zzz_zz_0_xz[i] * a_exp;

        g_z_0_0_0_zz_zz_0_yy[i] = -2.0 * g_z_zz_0_yy[i] + 2.0 * g_zzz_zz_0_yy[i] * a_exp;

        g_z_0_0_0_zz_zz_0_yz[i] = -2.0 * g_z_zz_0_yz[i] + 2.0 * g_zzz_zz_0_yz[i] * a_exp;

        g_z_0_0_0_zz_zz_0_zz[i] = -2.0 * g_z_zz_0_zz[i] + 2.0 * g_zzz_zz_0_zz[i] * a_exp;
    }
}

} // t4c_geom namespace

