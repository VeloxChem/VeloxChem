#include "GeomDeriv2000OfScalarForPPPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_pppd_0(CSimdArray<double>& buffer_2000_pppd,
                     const CSimdArray<double>& buffer_pppd,
                     const CSimdArray<double>& buffer_fppd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_pppd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pppd

    auto g_x_x_x_xx = buffer_pppd[0];

    auto g_x_x_x_xy = buffer_pppd[1];

    auto g_x_x_x_xz = buffer_pppd[2];

    auto g_x_x_x_yy = buffer_pppd[3];

    auto g_x_x_x_yz = buffer_pppd[4];

    auto g_x_x_x_zz = buffer_pppd[5];

    auto g_x_x_y_xx = buffer_pppd[6];

    auto g_x_x_y_xy = buffer_pppd[7];

    auto g_x_x_y_xz = buffer_pppd[8];

    auto g_x_x_y_yy = buffer_pppd[9];

    auto g_x_x_y_yz = buffer_pppd[10];

    auto g_x_x_y_zz = buffer_pppd[11];

    auto g_x_x_z_xx = buffer_pppd[12];

    auto g_x_x_z_xy = buffer_pppd[13];

    auto g_x_x_z_xz = buffer_pppd[14];

    auto g_x_x_z_yy = buffer_pppd[15];

    auto g_x_x_z_yz = buffer_pppd[16];

    auto g_x_x_z_zz = buffer_pppd[17];

    auto g_x_y_x_xx = buffer_pppd[18];

    auto g_x_y_x_xy = buffer_pppd[19];

    auto g_x_y_x_xz = buffer_pppd[20];

    auto g_x_y_x_yy = buffer_pppd[21];

    auto g_x_y_x_yz = buffer_pppd[22];

    auto g_x_y_x_zz = buffer_pppd[23];

    auto g_x_y_y_xx = buffer_pppd[24];

    auto g_x_y_y_xy = buffer_pppd[25];

    auto g_x_y_y_xz = buffer_pppd[26];

    auto g_x_y_y_yy = buffer_pppd[27];

    auto g_x_y_y_yz = buffer_pppd[28];

    auto g_x_y_y_zz = buffer_pppd[29];

    auto g_x_y_z_xx = buffer_pppd[30];

    auto g_x_y_z_xy = buffer_pppd[31];

    auto g_x_y_z_xz = buffer_pppd[32];

    auto g_x_y_z_yy = buffer_pppd[33];

    auto g_x_y_z_yz = buffer_pppd[34];

    auto g_x_y_z_zz = buffer_pppd[35];

    auto g_x_z_x_xx = buffer_pppd[36];

    auto g_x_z_x_xy = buffer_pppd[37];

    auto g_x_z_x_xz = buffer_pppd[38];

    auto g_x_z_x_yy = buffer_pppd[39];

    auto g_x_z_x_yz = buffer_pppd[40];

    auto g_x_z_x_zz = buffer_pppd[41];

    auto g_x_z_y_xx = buffer_pppd[42];

    auto g_x_z_y_xy = buffer_pppd[43];

    auto g_x_z_y_xz = buffer_pppd[44];

    auto g_x_z_y_yy = buffer_pppd[45];

    auto g_x_z_y_yz = buffer_pppd[46];

    auto g_x_z_y_zz = buffer_pppd[47];

    auto g_x_z_z_xx = buffer_pppd[48];

    auto g_x_z_z_xy = buffer_pppd[49];

    auto g_x_z_z_xz = buffer_pppd[50];

    auto g_x_z_z_yy = buffer_pppd[51];

    auto g_x_z_z_yz = buffer_pppd[52];

    auto g_x_z_z_zz = buffer_pppd[53];

    auto g_y_x_x_xx = buffer_pppd[54];

    auto g_y_x_x_xy = buffer_pppd[55];

    auto g_y_x_x_xz = buffer_pppd[56];

    auto g_y_x_x_yy = buffer_pppd[57];

    auto g_y_x_x_yz = buffer_pppd[58];

    auto g_y_x_x_zz = buffer_pppd[59];

    auto g_y_x_y_xx = buffer_pppd[60];

    auto g_y_x_y_xy = buffer_pppd[61];

    auto g_y_x_y_xz = buffer_pppd[62];

    auto g_y_x_y_yy = buffer_pppd[63];

    auto g_y_x_y_yz = buffer_pppd[64];

    auto g_y_x_y_zz = buffer_pppd[65];

    auto g_y_x_z_xx = buffer_pppd[66];

    auto g_y_x_z_xy = buffer_pppd[67];

    auto g_y_x_z_xz = buffer_pppd[68];

    auto g_y_x_z_yy = buffer_pppd[69];

    auto g_y_x_z_yz = buffer_pppd[70];

    auto g_y_x_z_zz = buffer_pppd[71];

    auto g_y_y_x_xx = buffer_pppd[72];

    auto g_y_y_x_xy = buffer_pppd[73];

    auto g_y_y_x_xz = buffer_pppd[74];

    auto g_y_y_x_yy = buffer_pppd[75];

    auto g_y_y_x_yz = buffer_pppd[76];

    auto g_y_y_x_zz = buffer_pppd[77];

    auto g_y_y_y_xx = buffer_pppd[78];

    auto g_y_y_y_xy = buffer_pppd[79];

    auto g_y_y_y_xz = buffer_pppd[80];

    auto g_y_y_y_yy = buffer_pppd[81];

    auto g_y_y_y_yz = buffer_pppd[82];

    auto g_y_y_y_zz = buffer_pppd[83];

    auto g_y_y_z_xx = buffer_pppd[84];

    auto g_y_y_z_xy = buffer_pppd[85];

    auto g_y_y_z_xz = buffer_pppd[86];

    auto g_y_y_z_yy = buffer_pppd[87];

    auto g_y_y_z_yz = buffer_pppd[88];

    auto g_y_y_z_zz = buffer_pppd[89];

    auto g_y_z_x_xx = buffer_pppd[90];

    auto g_y_z_x_xy = buffer_pppd[91];

    auto g_y_z_x_xz = buffer_pppd[92];

    auto g_y_z_x_yy = buffer_pppd[93];

    auto g_y_z_x_yz = buffer_pppd[94];

    auto g_y_z_x_zz = buffer_pppd[95];

    auto g_y_z_y_xx = buffer_pppd[96];

    auto g_y_z_y_xy = buffer_pppd[97];

    auto g_y_z_y_xz = buffer_pppd[98];

    auto g_y_z_y_yy = buffer_pppd[99];

    auto g_y_z_y_yz = buffer_pppd[100];

    auto g_y_z_y_zz = buffer_pppd[101];

    auto g_y_z_z_xx = buffer_pppd[102];

    auto g_y_z_z_xy = buffer_pppd[103];

    auto g_y_z_z_xz = buffer_pppd[104];

    auto g_y_z_z_yy = buffer_pppd[105];

    auto g_y_z_z_yz = buffer_pppd[106];

    auto g_y_z_z_zz = buffer_pppd[107];

    auto g_z_x_x_xx = buffer_pppd[108];

    auto g_z_x_x_xy = buffer_pppd[109];

    auto g_z_x_x_xz = buffer_pppd[110];

    auto g_z_x_x_yy = buffer_pppd[111];

    auto g_z_x_x_yz = buffer_pppd[112];

    auto g_z_x_x_zz = buffer_pppd[113];

    auto g_z_x_y_xx = buffer_pppd[114];

    auto g_z_x_y_xy = buffer_pppd[115];

    auto g_z_x_y_xz = buffer_pppd[116];

    auto g_z_x_y_yy = buffer_pppd[117];

    auto g_z_x_y_yz = buffer_pppd[118];

    auto g_z_x_y_zz = buffer_pppd[119];

    auto g_z_x_z_xx = buffer_pppd[120];

    auto g_z_x_z_xy = buffer_pppd[121];

    auto g_z_x_z_xz = buffer_pppd[122];

    auto g_z_x_z_yy = buffer_pppd[123];

    auto g_z_x_z_yz = buffer_pppd[124];

    auto g_z_x_z_zz = buffer_pppd[125];

    auto g_z_y_x_xx = buffer_pppd[126];

    auto g_z_y_x_xy = buffer_pppd[127];

    auto g_z_y_x_xz = buffer_pppd[128];

    auto g_z_y_x_yy = buffer_pppd[129];

    auto g_z_y_x_yz = buffer_pppd[130];

    auto g_z_y_x_zz = buffer_pppd[131];

    auto g_z_y_y_xx = buffer_pppd[132];

    auto g_z_y_y_xy = buffer_pppd[133];

    auto g_z_y_y_xz = buffer_pppd[134];

    auto g_z_y_y_yy = buffer_pppd[135];

    auto g_z_y_y_yz = buffer_pppd[136];

    auto g_z_y_y_zz = buffer_pppd[137];

    auto g_z_y_z_xx = buffer_pppd[138];

    auto g_z_y_z_xy = buffer_pppd[139];

    auto g_z_y_z_xz = buffer_pppd[140];

    auto g_z_y_z_yy = buffer_pppd[141];

    auto g_z_y_z_yz = buffer_pppd[142];

    auto g_z_y_z_zz = buffer_pppd[143];

    auto g_z_z_x_xx = buffer_pppd[144];

    auto g_z_z_x_xy = buffer_pppd[145];

    auto g_z_z_x_xz = buffer_pppd[146];

    auto g_z_z_x_yy = buffer_pppd[147];

    auto g_z_z_x_yz = buffer_pppd[148];

    auto g_z_z_x_zz = buffer_pppd[149];

    auto g_z_z_y_xx = buffer_pppd[150];

    auto g_z_z_y_xy = buffer_pppd[151];

    auto g_z_z_y_xz = buffer_pppd[152];

    auto g_z_z_y_yy = buffer_pppd[153];

    auto g_z_z_y_yz = buffer_pppd[154];

    auto g_z_z_y_zz = buffer_pppd[155];

    auto g_z_z_z_xx = buffer_pppd[156];

    auto g_z_z_z_xy = buffer_pppd[157];

    auto g_z_z_z_xz = buffer_pppd[158];

    auto g_z_z_z_yy = buffer_pppd[159];

    auto g_z_z_z_yz = buffer_pppd[160];

    auto g_z_z_z_zz = buffer_pppd[161];

    /// Set up components of auxilary buffer : buffer_fppd

    auto g_xxx_x_x_xx = buffer_fppd[0];

    auto g_xxx_x_x_xy = buffer_fppd[1];

    auto g_xxx_x_x_xz = buffer_fppd[2];

    auto g_xxx_x_x_yy = buffer_fppd[3];

    auto g_xxx_x_x_yz = buffer_fppd[4];

    auto g_xxx_x_x_zz = buffer_fppd[5];

    auto g_xxx_x_y_xx = buffer_fppd[6];

    auto g_xxx_x_y_xy = buffer_fppd[7];

    auto g_xxx_x_y_xz = buffer_fppd[8];

    auto g_xxx_x_y_yy = buffer_fppd[9];

    auto g_xxx_x_y_yz = buffer_fppd[10];

    auto g_xxx_x_y_zz = buffer_fppd[11];

    auto g_xxx_x_z_xx = buffer_fppd[12];

    auto g_xxx_x_z_xy = buffer_fppd[13];

    auto g_xxx_x_z_xz = buffer_fppd[14];

    auto g_xxx_x_z_yy = buffer_fppd[15];

    auto g_xxx_x_z_yz = buffer_fppd[16];

    auto g_xxx_x_z_zz = buffer_fppd[17];

    auto g_xxx_y_x_xx = buffer_fppd[18];

    auto g_xxx_y_x_xy = buffer_fppd[19];

    auto g_xxx_y_x_xz = buffer_fppd[20];

    auto g_xxx_y_x_yy = buffer_fppd[21];

    auto g_xxx_y_x_yz = buffer_fppd[22];

    auto g_xxx_y_x_zz = buffer_fppd[23];

    auto g_xxx_y_y_xx = buffer_fppd[24];

    auto g_xxx_y_y_xy = buffer_fppd[25];

    auto g_xxx_y_y_xz = buffer_fppd[26];

    auto g_xxx_y_y_yy = buffer_fppd[27];

    auto g_xxx_y_y_yz = buffer_fppd[28];

    auto g_xxx_y_y_zz = buffer_fppd[29];

    auto g_xxx_y_z_xx = buffer_fppd[30];

    auto g_xxx_y_z_xy = buffer_fppd[31];

    auto g_xxx_y_z_xz = buffer_fppd[32];

    auto g_xxx_y_z_yy = buffer_fppd[33];

    auto g_xxx_y_z_yz = buffer_fppd[34];

    auto g_xxx_y_z_zz = buffer_fppd[35];

    auto g_xxx_z_x_xx = buffer_fppd[36];

    auto g_xxx_z_x_xy = buffer_fppd[37];

    auto g_xxx_z_x_xz = buffer_fppd[38];

    auto g_xxx_z_x_yy = buffer_fppd[39];

    auto g_xxx_z_x_yz = buffer_fppd[40];

    auto g_xxx_z_x_zz = buffer_fppd[41];

    auto g_xxx_z_y_xx = buffer_fppd[42];

    auto g_xxx_z_y_xy = buffer_fppd[43];

    auto g_xxx_z_y_xz = buffer_fppd[44];

    auto g_xxx_z_y_yy = buffer_fppd[45];

    auto g_xxx_z_y_yz = buffer_fppd[46];

    auto g_xxx_z_y_zz = buffer_fppd[47];

    auto g_xxx_z_z_xx = buffer_fppd[48];

    auto g_xxx_z_z_xy = buffer_fppd[49];

    auto g_xxx_z_z_xz = buffer_fppd[50];

    auto g_xxx_z_z_yy = buffer_fppd[51];

    auto g_xxx_z_z_yz = buffer_fppd[52];

    auto g_xxx_z_z_zz = buffer_fppd[53];

    auto g_xxy_x_x_xx = buffer_fppd[54];

    auto g_xxy_x_x_xy = buffer_fppd[55];

    auto g_xxy_x_x_xz = buffer_fppd[56];

    auto g_xxy_x_x_yy = buffer_fppd[57];

    auto g_xxy_x_x_yz = buffer_fppd[58];

    auto g_xxy_x_x_zz = buffer_fppd[59];

    auto g_xxy_x_y_xx = buffer_fppd[60];

    auto g_xxy_x_y_xy = buffer_fppd[61];

    auto g_xxy_x_y_xz = buffer_fppd[62];

    auto g_xxy_x_y_yy = buffer_fppd[63];

    auto g_xxy_x_y_yz = buffer_fppd[64];

    auto g_xxy_x_y_zz = buffer_fppd[65];

    auto g_xxy_x_z_xx = buffer_fppd[66];

    auto g_xxy_x_z_xy = buffer_fppd[67];

    auto g_xxy_x_z_xz = buffer_fppd[68];

    auto g_xxy_x_z_yy = buffer_fppd[69];

    auto g_xxy_x_z_yz = buffer_fppd[70];

    auto g_xxy_x_z_zz = buffer_fppd[71];

    auto g_xxy_y_x_xx = buffer_fppd[72];

    auto g_xxy_y_x_xy = buffer_fppd[73];

    auto g_xxy_y_x_xz = buffer_fppd[74];

    auto g_xxy_y_x_yy = buffer_fppd[75];

    auto g_xxy_y_x_yz = buffer_fppd[76];

    auto g_xxy_y_x_zz = buffer_fppd[77];

    auto g_xxy_y_y_xx = buffer_fppd[78];

    auto g_xxy_y_y_xy = buffer_fppd[79];

    auto g_xxy_y_y_xz = buffer_fppd[80];

    auto g_xxy_y_y_yy = buffer_fppd[81];

    auto g_xxy_y_y_yz = buffer_fppd[82];

    auto g_xxy_y_y_zz = buffer_fppd[83];

    auto g_xxy_y_z_xx = buffer_fppd[84];

    auto g_xxy_y_z_xy = buffer_fppd[85];

    auto g_xxy_y_z_xz = buffer_fppd[86];

    auto g_xxy_y_z_yy = buffer_fppd[87];

    auto g_xxy_y_z_yz = buffer_fppd[88];

    auto g_xxy_y_z_zz = buffer_fppd[89];

    auto g_xxy_z_x_xx = buffer_fppd[90];

    auto g_xxy_z_x_xy = buffer_fppd[91];

    auto g_xxy_z_x_xz = buffer_fppd[92];

    auto g_xxy_z_x_yy = buffer_fppd[93];

    auto g_xxy_z_x_yz = buffer_fppd[94];

    auto g_xxy_z_x_zz = buffer_fppd[95];

    auto g_xxy_z_y_xx = buffer_fppd[96];

    auto g_xxy_z_y_xy = buffer_fppd[97];

    auto g_xxy_z_y_xz = buffer_fppd[98];

    auto g_xxy_z_y_yy = buffer_fppd[99];

    auto g_xxy_z_y_yz = buffer_fppd[100];

    auto g_xxy_z_y_zz = buffer_fppd[101];

    auto g_xxy_z_z_xx = buffer_fppd[102];

    auto g_xxy_z_z_xy = buffer_fppd[103];

    auto g_xxy_z_z_xz = buffer_fppd[104];

    auto g_xxy_z_z_yy = buffer_fppd[105];

    auto g_xxy_z_z_yz = buffer_fppd[106];

    auto g_xxy_z_z_zz = buffer_fppd[107];

    auto g_xxz_x_x_xx = buffer_fppd[108];

    auto g_xxz_x_x_xy = buffer_fppd[109];

    auto g_xxz_x_x_xz = buffer_fppd[110];

    auto g_xxz_x_x_yy = buffer_fppd[111];

    auto g_xxz_x_x_yz = buffer_fppd[112];

    auto g_xxz_x_x_zz = buffer_fppd[113];

    auto g_xxz_x_y_xx = buffer_fppd[114];

    auto g_xxz_x_y_xy = buffer_fppd[115];

    auto g_xxz_x_y_xz = buffer_fppd[116];

    auto g_xxz_x_y_yy = buffer_fppd[117];

    auto g_xxz_x_y_yz = buffer_fppd[118];

    auto g_xxz_x_y_zz = buffer_fppd[119];

    auto g_xxz_x_z_xx = buffer_fppd[120];

    auto g_xxz_x_z_xy = buffer_fppd[121];

    auto g_xxz_x_z_xz = buffer_fppd[122];

    auto g_xxz_x_z_yy = buffer_fppd[123];

    auto g_xxz_x_z_yz = buffer_fppd[124];

    auto g_xxz_x_z_zz = buffer_fppd[125];

    auto g_xxz_y_x_xx = buffer_fppd[126];

    auto g_xxz_y_x_xy = buffer_fppd[127];

    auto g_xxz_y_x_xz = buffer_fppd[128];

    auto g_xxz_y_x_yy = buffer_fppd[129];

    auto g_xxz_y_x_yz = buffer_fppd[130];

    auto g_xxz_y_x_zz = buffer_fppd[131];

    auto g_xxz_y_y_xx = buffer_fppd[132];

    auto g_xxz_y_y_xy = buffer_fppd[133];

    auto g_xxz_y_y_xz = buffer_fppd[134];

    auto g_xxz_y_y_yy = buffer_fppd[135];

    auto g_xxz_y_y_yz = buffer_fppd[136];

    auto g_xxz_y_y_zz = buffer_fppd[137];

    auto g_xxz_y_z_xx = buffer_fppd[138];

    auto g_xxz_y_z_xy = buffer_fppd[139];

    auto g_xxz_y_z_xz = buffer_fppd[140];

    auto g_xxz_y_z_yy = buffer_fppd[141];

    auto g_xxz_y_z_yz = buffer_fppd[142];

    auto g_xxz_y_z_zz = buffer_fppd[143];

    auto g_xxz_z_x_xx = buffer_fppd[144];

    auto g_xxz_z_x_xy = buffer_fppd[145];

    auto g_xxz_z_x_xz = buffer_fppd[146];

    auto g_xxz_z_x_yy = buffer_fppd[147];

    auto g_xxz_z_x_yz = buffer_fppd[148];

    auto g_xxz_z_x_zz = buffer_fppd[149];

    auto g_xxz_z_y_xx = buffer_fppd[150];

    auto g_xxz_z_y_xy = buffer_fppd[151];

    auto g_xxz_z_y_xz = buffer_fppd[152];

    auto g_xxz_z_y_yy = buffer_fppd[153];

    auto g_xxz_z_y_yz = buffer_fppd[154];

    auto g_xxz_z_y_zz = buffer_fppd[155];

    auto g_xxz_z_z_xx = buffer_fppd[156];

    auto g_xxz_z_z_xy = buffer_fppd[157];

    auto g_xxz_z_z_xz = buffer_fppd[158];

    auto g_xxz_z_z_yy = buffer_fppd[159];

    auto g_xxz_z_z_yz = buffer_fppd[160];

    auto g_xxz_z_z_zz = buffer_fppd[161];

    auto g_xyy_x_x_xx = buffer_fppd[162];

    auto g_xyy_x_x_xy = buffer_fppd[163];

    auto g_xyy_x_x_xz = buffer_fppd[164];

    auto g_xyy_x_x_yy = buffer_fppd[165];

    auto g_xyy_x_x_yz = buffer_fppd[166];

    auto g_xyy_x_x_zz = buffer_fppd[167];

    auto g_xyy_x_y_xx = buffer_fppd[168];

    auto g_xyy_x_y_xy = buffer_fppd[169];

    auto g_xyy_x_y_xz = buffer_fppd[170];

    auto g_xyy_x_y_yy = buffer_fppd[171];

    auto g_xyy_x_y_yz = buffer_fppd[172];

    auto g_xyy_x_y_zz = buffer_fppd[173];

    auto g_xyy_x_z_xx = buffer_fppd[174];

    auto g_xyy_x_z_xy = buffer_fppd[175];

    auto g_xyy_x_z_xz = buffer_fppd[176];

    auto g_xyy_x_z_yy = buffer_fppd[177];

    auto g_xyy_x_z_yz = buffer_fppd[178];

    auto g_xyy_x_z_zz = buffer_fppd[179];

    auto g_xyy_y_x_xx = buffer_fppd[180];

    auto g_xyy_y_x_xy = buffer_fppd[181];

    auto g_xyy_y_x_xz = buffer_fppd[182];

    auto g_xyy_y_x_yy = buffer_fppd[183];

    auto g_xyy_y_x_yz = buffer_fppd[184];

    auto g_xyy_y_x_zz = buffer_fppd[185];

    auto g_xyy_y_y_xx = buffer_fppd[186];

    auto g_xyy_y_y_xy = buffer_fppd[187];

    auto g_xyy_y_y_xz = buffer_fppd[188];

    auto g_xyy_y_y_yy = buffer_fppd[189];

    auto g_xyy_y_y_yz = buffer_fppd[190];

    auto g_xyy_y_y_zz = buffer_fppd[191];

    auto g_xyy_y_z_xx = buffer_fppd[192];

    auto g_xyy_y_z_xy = buffer_fppd[193];

    auto g_xyy_y_z_xz = buffer_fppd[194];

    auto g_xyy_y_z_yy = buffer_fppd[195];

    auto g_xyy_y_z_yz = buffer_fppd[196];

    auto g_xyy_y_z_zz = buffer_fppd[197];

    auto g_xyy_z_x_xx = buffer_fppd[198];

    auto g_xyy_z_x_xy = buffer_fppd[199];

    auto g_xyy_z_x_xz = buffer_fppd[200];

    auto g_xyy_z_x_yy = buffer_fppd[201];

    auto g_xyy_z_x_yz = buffer_fppd[202];

    auto g_xyy_z_x_zz = buffer_fppd[203];

    auto g_xyy_z_y_xx = buffer_fppd[204];

    auto g_xyy_z_y_xy = buffer_fppd[205];

    auto g_xyy_z_y_xz = buffer_fppd[206];

    auto g_xyy_z_y_yy = buffer_fppd[207];

    auto g_xyy_z_y_yz = buffer_fppd[208];

    auto g_xyy_z_y_zz = buffer_fppd[209];

    auto g_xyy_z_z_xx = buffer_fppd[210];

    auto g_xyy_z_z_xy = buffer_fppd[211];

    auto g_xyy_z_z_xz = buffer_fppd[212];

    auto g_xyy_z_z_yy = buffer_fppd[213];

    auto g_xyy_z_z_yz = buffer_fppd[214];

    auto g_xyy_z_z_zz = buffer_fppd[215];

    auto g_xyz_x_x_xx = buffer_fppd[216];

    auto g_xyz_x_x_xy = buffer_fppd[217];

    auto g_xyz_x_x_xz = buffer_fppd[218];

    auto g_xyz_x_x_yy = buffer_fppd[219];

    auto g_xyz_x_x_yz = buffer_fppd[220];

    auto g_xyz_x_x_zz = buffer_fppd[221];

    auto g_xyz_x_y_xx = buffer_fppd[222];

    auto g_xyz_x_y_xy = buffer_fppd[223];

    auto g_xyz_x_y_xz = buffer_fppd[224];

    auto g_xyz_x_y_yy = buffer_fppd[225];

    auto g_xyz_x_y_yz = buffer_fppd[226];

    auto g_xyz_x_y_zz = buffer_fppd[227];

    auto g_xyz_x_z_xx = buffer_fppd[228];

    auto g_xyz_x_z_xy = buffer_fppd[229];

    auto g_xyz_x_z_xz = buffer_fppd[230];

    auto g_xyz_x_z_yy = buffer_fppd[231];

    auto g_xyz_x_z_yz = buffer_fppd[232];

    auto g_xyz_x_z_zz = buffer_fppd[233];

    auto g_xyz_y_x_xx = buffer_fppd[234];

    auto g_xyz_y_x_xy = buffer_fppd[235];

    auto g_xyz_y_x_xz = buffer_fppd[236];

    auto g_xyz_y_x_yy = buffer_fppd[237];

    auto g_xyz_y_x_yz = buffer_fppd[238];

    auto g_xyz_y_x_zz = buffer_fppd[239];

    auto g_xyz_y_y_xx = buffer_fppd[240];

    auto g_xyz_y_y_xy = buffer_fppd[241];

    auto g_xyz_y_y_xz = buffer_fppd[242];

    auto g_xyz_y_y_yy = buffer_fppd[243];

    auto g_xyz_y_y_yz = buffer_fppd[244];

    auto g_xyz_y_y_zz = buffer_fppd[245];

    auto g_xyz_y_z_xx = buffer_fppd[246];

    auto g_xyz_y_z_xy = buffer_fppd[247];

    auto g_xyz_y_z_xz = buffer_fppd[248];

    auto g_xyz_y_z_yy = buffer_fppd[249];

    auto g_xyz_y_z_yz = buffer_fppd[250];

    auto g_xyz_y_z_zz = buffer_fppd[251];

    auto g_xyz_z_x_xx = buffer_fppd[252];

    auto g_xyz_z_x_xy = buffer_fppd[253];

    auto g_xyz_z_x_xz = buffer_fppd[254];

    auto g_xyz_z_x_yy = buffer_fppd[255];

    auto g_xyz_z_x_yz = buffer_fppd[256];

    auto g_xyz_z_x_zz = buffer_fppd[257];

    auto g_xyz_z_y_xx = buffer_fppd[258];

    auto g_xyz_z_y_xy = buffer_fppd[259];

    auto g_xyz_z_y_xz = buffer_fppd[260];

    auto g_xyz_z_y_yy = buffer_fppd[261];

    auto g_xyz_z_y_yz = buffer_fppd[262];

    auto g_xyz_z_y_zz = buffer_fppd[263];

    auto g_xyz_z_z_xx = buffer_fppd[264];

    auto g_xyz_z_z_xy = buffer_fppd[265];

    auto g_xyz_z_z_xz = buffer_fppd[266];

    auto g_xyz_z_z_yy = buffer_fppd[267];

    auto g_xyz_z_z_yz = buffer_fppd[268];

    auto g_xyz_z_z_zz = buffer_fppd[269];

    auto g_xzz_x_x_xx = buffer_fppd[270];

    auto g_xzz_x_x_xy = buffer_fppd[271];

    auto g_xzz_x_x_xz = buffer_fppd[272];

    auto g_xzz_x_x_yy = buffer_fppd[273];

    auto g_xzz_x_x_yz = buffer_fppd[274];

    auto g_xzz_x_x_zz = buffer_fppd[275];

    auto g_xzz_x_y_xx = buffer_fppd[276];

    auto g_xzz_x_y_xy = buffer_fppd[277];

    auto g_xzz_x_y_xz = buffer_fppd[278];

    auto g_xzz_x_y_yy = buffer_fppd[279];

    auto g_xzz_x_y_yz = buffer_fppd[280];

    auto g_xzz_x_y_zz = buffer_fppd[281];

    auto g_xzz_x_z_xx = buffer_fppd[282];

    auto g_xzz_x_z_xy = buffer_fppd[283];

    auto g_xzz_x_z_xz = buffer_fppd[284];

    auto g_xzz_x_z_yy = buffer_fppd[285];

    auto g_xzz_x_z_yz = buffer_fppd[286];

    auto g_xzz_x_z_zz = buffer_fppd[287];

    auto g_xzz_y_x_xx = buffer_fppd[288];

    auto g_xzz_y_x_xy = buffer_fppd[289];

    auto g_xzz_y_x_xz = buffer_fppd[290];

    auto g_xzz_y_x_yy = buffer_fppd[291];

    auto g_xzz_y_x_yz = buffer_fppd[292];

    auto g_xzz_y_x_zz = buffer_fppd[293];

    auto g_xzz_y_y_xx = buffer_fppd[294];

    auto g_xzz_y_y_xy = buffer_fppd[295];

    auto g_xzz_y_y_xz = buffer_fppd[296];

    auto g_xzz_y_y_yy = buffer_fppd[297];

    auto g_xzz_y_y_yz = buffer_fppd[298];

    auto g_xzz_y_y_zz = buffer_fppd[299];

    auto g_xzz_y_z_xx = buffer_fppd[300];

    auto g_xzz_y_z_xy = buffer_fppd[301];

    auto g_xzz_y_z_xz = buffer_fppd[302];

    auto g_xzz_y_z_yy = buffer_fppd[303];

    auto g_xzz_y_z_yz = buffer_fppd[304];

    auto g_xzz_y_z_zz = buffer_fppd[305];

    auto g_xzz_z_x_xx = buffer_fppd[306];

    auto g_xzz_z_x_xy = buffer_fppd[307];

    auto g_xzz_z_x_xz = buffer_fppd[308];

    auto g_xzz_z_x_yy = buffer_fppd[309];

    auto g_xzz_z_x_yz = buffer_fppd[310];

    auto g_xzz_z_x_zz = buffer_fppd[311];

    auto g_xzz_z_y_xx = buffer_fppd[312];

    auto g_xzz_z_y_xy = buffer_fppd[313];

    auto g_xzz_z_y_xz = buffer_fppd[314];

    auto g_xzz_z_y_yy = buffer_fppd[315];

    auto g_xzz_z_y_yz = buffer_fppd[316];

    auto g_xzz_z_y_zz = buffer_fppd[317];

    auto g_xzz_z_z_xx = buffer_fppd[318];

    auto g_xzz_z_z_xy = buffer_fppd[319];

    auto g_xzz_z_z_xz = buffer_fppd[320];

    auto g_xzz_z_z_yy = buffer_fppd[321];

    auto g_xzz_z_z_yz = buffer_fppd[322];

    auto g_xzz_z_z_zz = buffer_fppd[323];

    auto g_yyy_x_x_xx = buffer_fppd[324];

    auto g_yyy_x_x_xy = buffer_fppd[325];

    auto g_yyy_x_x_xz = buffer_fppd[326];

    auto g_yyy_x_x_yy = buffer_fppd[327];

    auto g_yyy_x_x_yz = buffer_fppd[328];

    auto g_yyy_x_x_zz = buffer_fppd[329];

    auto g_yyy_x_y_xx = buffer_fppd[330];

    auto g_yyy_x_y_xy = buffer_fppd[331];

    auto g_yyy_x_y_xz = buffer_fppd[332];

    auto g_yyy_x_y_yy = buffer_fppd[333];

    auto g_yyy_x_y_yz = buffer_fppd[334];

    auto g_yyy_x_y_zz = buffer_fppd[335];

    auto g_yyy_x_z_xx = buffer_fppd[336];

    auto g_yyy_x_z_xy = buffer_fppd[337];

    auto g_yyy_x_z_xz = buffer_fppd[338];

    auto g_yyy_x_z_yy = buffer_fppd[339];

    auto g_yyy_x_z_yz = buffer_fppd[340];

    auto g_yyy_x_z_zz = buffer_fppd[341];

    auto g_yyy_y_x_xx = buffer_fppd[342];

    auto g_yyy_y_x_xy = buffer_fppd[343];

    auto g_yyy_y_x_xz = buffer_fppd[344];

    auto g_yyy_y_x_yy = buffer_fppd[345];

    auto g_yyy_y_x_yz = buffer_fppd[346];

    auto g_yyy_y_x_zz = buffer_fppd[347];

    auto g_yyy_y_y_xx = buffer_fppd[348];

    auto g_yyy_y_y_xy = buffer_fppd[349];

    auto g_yyy_y_y_xz = buffer_fppd[350];

    auto g_yyy_y_y_yy = buffer_fppd[351];

    auto g_yyy_y_y_yz = buffer_fppd[352];

    auto g_yyy_y_y_zz = buffer_fppd[353];

    auto g_yyy_y_z_xx = buffer_fppd[354];

    auto g_yyy_y_z_xy = buffer_fppd[355];

    auto g_yyy_y_z_xz = buffer_fppd[356];

    auto g_yyy_y_z_yy = buffer_fppd[357];

    auto g_yyy_y_z_yz = buffer_fppd[358];

    auto g_yyy_y_z_zz = buffer_fppd[359];

    auto g_yyy_z_x_xx = buffer_fppd[360];

    auto g_yyy_z_x_xy = buffer_fppd[361];

    auto g_yyy_z_x_xz = buffer_fppd[362];

    auto g_yyy_z_x_yy = buffer_fppd[363];

    auto g_yyy_z_x_yz = buffer_fppd[364];

    auto g_yyy_z_x_zz = buffer_fppd[365];

    auto g_yyy_z_y_xx = buffer_fppd[366];

    auto g_yyy_z_y_xy = buffer_fppd[367];

    auto g_yyy_z_y_xz = buffer_fppd[368];

    auto g_yyy_z_y_yy = buffer_fppd[369];

    auto g_yyy_z_y_yz = buffer_fppd[370];

    auto g_yyy_z_y_zz = buffer_fppd[371];

    auto g_yyy_z_z_xx = buffer_fppd[372];

    auto g_yyy_z_z_xy = buffer_fppd[373];

    auto g_yyy_z_z_xz = buffer_fppd[374];

    auto g_yyy_z_z_yy = buffer_fppd[375];

    auto g_yyy_z_z_yz = buffer_fppd[376];

    auto g_yyy_z_z_zz = buffer_fppd[377];

    auto g_yyz_x_x_xx = buffer_fppd[378];

    auto g_yyz_x_x_xy = buffer_fppd[379];

    auto g_yyz_x_x_xz = buffer_fppd[380];

    auto g_yyz_x_x_yy = buffer_fppd[381];

    auto g_yyz_x_x_yz = buffer_fppd[382];

    auto g_yyz_x_x_zz = buffer_fppd[383];

    auto g_yyz_x_y_xx = buffer_fppd[384];

    auto g_yyz_x_y_xy = buffer_fppd[385];

    auto g_yyz_x_y_xz = buffer_fppd[386];

    auto g_yyz_x_y_yy = buffer_fppd[387];

    auto g_yyz_x_y_yz = buffer_fppd[388];

    auto g_yyz_x_y_zz = buffer_fppd[389];

    auto g_yyz_x_z_xx = buffer_fppd[390];

    auto g_yyz_x_z_xy = buffer_fppd[391];

    auto g_yyz_x_z_xz = buffer_fppd[392];

    auto g_yyz_x_z_yy = buffer_fppd[393];

    auto g_yyz_x_z_yz = buffer_fppd[394];

    auto g_yyz_x_z_zz = buffer_fppd[395];

    auto g_yyz_y_x_xx = buffer_fppd[396];

    auto g_yyz_y_x_xy = buffer_fppd[397];

    auto g_yyz_y_x_xz = buffer_fppd[398];

    auto g_yyz_y_x_yy = buffer_fppd[399];

    auto g_yyz_y_x_yz = buffer_fppd[400];

    auto g_yyz_y_x_zz = buffer_fppd[401];

    auto g_yyz_y_y_xx = buffer_fppd[402];

    auto g_yyz_y_y_xy = buffer_fppd[403];

    auto g_yyz_y_y_xz = buffer_fppd[404];

    auto g_yyz_y_y_yy = buffer_fppd[405];

    auto g_yyz_y_y_yz = buffer_fppd[406];

    auto g_yyz_y_y_zz = buffer_fppd[407];

    auto g_yyz_y_z_xx = buffer_fppd[408];

    auto g_yyz_y_z_xy = buffer_fppd[409];

    auto g_yyz_y_z_xz = buffer_fppd[410];

    auto g_yyz_y_z_yy = buffer_fppd[411];

    auto g_yyz_y_z_yz = buffer_fppd[412];

    auto g_yyz_y_z_zz = buffer_fppd[413];

    auto g_yyz_z_x_xx = buffer_fppd[414];

    auto g_yyz_z_x_xy = buffer_fppd[415];

    auto g_yyz_z_x_xz = buffer_fppd[416];

    auto g_yyz_z_x_yy = buffer_fppd[417];

    auto g_yyz_z_x_yz = buffer_fppd[418];

    auto g_yyz_z_x_zz = buffer_fppd[419];

    auto g_yyz_z_y_xx = buffer_fppd[420];

    auto g_yyz_z_y_xy = buffer_fppd[421];

    auto g_yyz_z_y_xz = buffer_fppd[422];

    auto g_yyz_z_y_yy = buffer_fppd[423];

    auto g_yyz_z_y_yz = buffer_fppd[424];

    auto g_yyz_z_y_zz = buffer_fppd[425];

    auto g_yyz_z_z_xx = buffer_fppd[426];

    auto g_yyz_z_z_xy = buffer_fppd[427];

    auto g_yyz_z_z_xz = buffer_fppd[428];

    auto g_yyz_z_z_yy = buffer_fppd[429];

    auto g_yyz_z_z_yz = buffer_fppd[430];

    auto g_yyz_z_z_zz = buffer_fppd[431];

    auto g_yzz_x_x_xx = buffer_fppd[432];

    auto g_yzz_x_x_xy = buffer_fppd[433];

    auto g_yzz_x_x_xz = buffer_fppd[434];

    auto g_yzz_x_x_yy = buffer_fppd[435];

    auto g_yzz_x_x_yz = buffer_fppd[436];

    auto g_yzz_x_x_zz = buffer_fppd[437];

    auto g_yzz_x_y_xx = buffer_fppd[438];

    auto g_yzz_x_y_xy = buffer_fppd[439];

    auto g_yzz_x_y_xz = buffer_fppd[440];

    auto g_yzz_x_y_yy = buffer_fppd[441];

    auto g_yzz_x_y_yz = buffer_fppd[442];

    auto g_yzz_x_y_zz = buffer_fppd[443];

    auto g_yzz_x_z_xx = buffer_fppd[444];

    auto g_yzz_x_z_xy = buffer_fppd[445];

    auto g_yzz_x_z_xz = buffer_fppd[446];

    auto g_yzz_x_z_yy = buffer_fppd[447];

    auto g_yzz_x_z_yz = buffer_fppd[448];

    auto g_yzz_x_z_zz = buffer_fppd[449];

    auto g_yzz_y_x_xx = buffer_fppd[450];

    auto g_yzz_y_x_xy = buffer_fppd[451];

    auto g_yzz_y_x_xz = buffer_fppd[452];

    auto g_yzz_y_x_yy = buffer_fppd[453];

    auto g_yzz_y_x_yz = buffer_fppd[454];

    auto g_yzz_y_x_zz = buffer_fppd[455];

    auto g_yzz_y_y_xx = buffer_fppd[456];

    auto g_yzz_y_y_xy = buffer_fppd[457];

    auto g_yzz_y_y_xz = buffer_fppd[458];

    auto g_yzz_y_y_yy = buffer_fppd[459];

    auto g_yzz_y_y_yz = buffer_fppd[460];

    auto g_yzz_y_y_zz = buffer_fppd[461];

    auto g_yzz_y_z_xx = buffer_fppd[462];

    auto g_yzz_y_z_xy = buffer_fppd[463];

    auto g_yzz_y_z_xz = buffer_fppd[464];

    auto g_yzz_y_z_yy = buffer_fppd[465];

    auto g_yzz_y_z_yz = buffer_fppd[466];

    auto g_yzz_y_z_zz = buffer_fppd[467];

    auto g_yzz_z_x_xx = buffer_fppd[468];

    auto g_yzz_z_x_xy = buffer_fppd[469];

    auto g_yzz_z_x_xz = buffer_fppd[470];

    auto g_yzz_z_x_yy = buffer_fppd[471];

    auto g_yzz_z_x_yz = buffer_fppd[472];

    auto g_yzz_z_x_zz = buffer_fppd[473];

    auto g_yzz_z_y_xx = buffer_fppd[474];

    auto g_yzz_z_y_xy = buffer_fppd[475];

    auto g_yzz_z_y_xz = buffer_fppd[476];

    auto g_yzz_z_y_yy = buffer_fppd[477];

    auto g_yzz_z_y_yz = buffer_fppd[478];

    auto g_yzz_z_y_zz = buffer_fppd[479];

    auto g_yzz_z_z_xx = buffer_fppd[480];

    auto g_yzz_z_z_xy = buffer_fppd[481];

    auto g_yzz_z_z_xz = buffer_fppd[482];

    auto g_yzz_z_z_yy = buffer_fppd[483];

    auto g_yzz_z_z_yz = buffer_fppd[484];

    auto g_yzz_z_z_zz = buffer_fppd[485];

    auto g_zzz_x_x_xx = buffer_fppd[486];

    auto g_zzz_x_x_xy = buffer_fppd[487];

    auto g_zzz_x_x_xz = buffer_fppd[488];

    auto g_zzz_x_x_yy = buffer_fppd[489];

    auto g_zzz_x_x_yz = buffer_fppd[490];

    auto g_zzz_x_x_zz = buffer_fppd[491];

    auto g_zzz_x_y_xx = buffer_fppd[492];

    auto g_zzz_x_y_xy = buffer_fppd[493];

    auto g_zzz_x_y_xz = buffer_fppd[494];

    auto g_zzz_x_y_yy = buffer_fppd[495];

    auto g_zzz_x_y_yz = buffer_fppd[496];

    auto g_zzz_x_y_zz = buffer_fppd[497];

    auto g_zzz_x_z_xx = buffer_fppd[498];

    auto g_zzz_x_z_xy = buffer_fppd[499];

    auto g_zzz_x_z_xz = buffer_fppd[500];

    auto g_zzz_x_z_yy = buffer_fppd[501];

    auto g_zzz_x_z_yz = buffer_fppd[502];

    auto g_zzz_x_z_zz = buffer_fppd[503];

    auto g_zzz_y_x_xx = buffer_fppd[504];

    auto g_zzz_y_x_xy = buffer_fppd[505];

    auto g_zzz_y_x_xz = buffer_fppd[506];

    auto g_zzz_y_x_yy = buffer_fppd[507];

    auto g_zzz_y_x_yz = buffer_fppd[508];

    auto g_zzz_y_x_zz = buffer_fppd[509];

    auto g_zzz_y_y_xx = buffer_fppd[510];

    auto g_zzz_y_y_xy = buffer_fppd[511];

    auto g_zzz_y_y_xz = buffer_fppd[512];

    auto g_zzz_y_y_yy = buffer_fppd[513];

    auto g_zzz_y_y_yz = buffer_fppd[514];

    auto g_zzz_y_y_zz = buffer_fppd[515];

    auto g_zzz_y_z_xx = buffer_fppd[516];

    auto g_zzz_y_z_xy = buffer_fppd[517];

    auto g_zzz_y_z_xz = buffer_fppd[518];

    auto g_zzz_y_z_yy = buffer_fppd[519];

    auto g_zzz_y_z_yz = buffer_fppd[520];

    auto g_zzz_y_z_zz = buffer_fppd[521];

    auto g_zzz_z_x_xx = buffer_fppd[522];

    auto g_zzz_z_x_xy = buffer_fppd[523];

    auto g_zzz_z_x_xz = buffer_fppd[524];

    auto g_zzz_z_x_yy = buffer_fppd[525];

    auto g_zzz_z_x_yz = buffer_fppd[526];

    auto g_zzz_z_x_zz = buffer_fppd[527];

    auto g_zzz_z_y_xx = buffer_fppd[528];

    auto g_zzz_z_y_xy = buffer_fppd[529];

    auto g_zzz_z_y_xz = buffer_fppd[530];

    auto g_zzz_z_y_yy = buffer_fppd[531];

    auto g_zzz_z_y_yz = buffer_fppd[532];

    auto g_zzz_z_y_zz = buffer_fppd[533];

    auto g_zzz_z_z_xx = buffer_fppd[534];

    auto g_zzz_z_z_xy = buffer_fppd[535];

    auto g_zzz_z_z_xz = buffer_fppd[536];

    auto g_zzz_z_z_yy = buffer_fppd[537];

    auto g_zzz_z_z_yz = buffer_fppd[538];

    auto g_zzz_z_z_zz = buffer_fppd[539];

    /// Set up components of integrals buffer : buffer_2000_pppd

    auto g_xx_0_0_0_x_x_x_xx = buffer_2000_pppd[0];

    auto g_xx_0_0_0_x_x_x_xy = buffer_2000_pppd[1];

    auto g_xx_0_0_0_x_x_x_xz = buffer_2000_pppd[2];

    auto g_xx_0_0_0_x_x_x_yy = buffer_2000_pppd[3];

    auto g_xx_0_0_0_x_x_x_yz = buffer_2000_pppd[4];

    auto g_xx_0_0_0_x_x_x_zz = buffer_2000_pppd[5];

    auto g_xx_0_0_0_x_x_y_xx = buffer_2000_pppd[6];

    auto g_xx_0_0_0_x_x_y_xy = buffer_2000_pppd[7];

    auto g_xx_0_0_0_x_x_y_xz = buffer_2000_pppd[8];

    auto g_xx_0_0_0_x_x_y_yy = buffer_2000_pppd[9];

    auto g_xx_0_0_0_x_x_y_yz = buffer_2000_pppd[10];

    auto g_xx_0_0_0_x_x_y_zz = buffer_2000_pppd[11];

    auto g_xx_0_0_0_x_x_z_xx = buffer_2000_pppd[12];

    auto g_xx_0_0_0_x_x_z_xy = buffer_2000_pppd[13];

    auto g_xx_0_0_0_x_x_z_xz = buffer_2000_pppd[14];

    auto g_xx_0_0_0_x_x_z_yy = buffer_2000_pppd[15];

    auto g_xx_0_0_0_x_x_z_yz = buffer_2000_pppd[16];

    auto g_xx_0_0_0_x_x_z_zz = buffer_2000_pppd[17];

    auto g_xx_0_0_0_x_y_x_xx = buffer_2000_pppd[18];

    auto g_xx_0_0_0_x_y_x_xy = buffer_2000_pppd[19];

    auto g_xx_0_0_0_x_y_x_xz = buffer_2000_pppd[20];

    auto g_xx_0_0_0_x_y_x_yy = buffer_2000_pppd[21];

    auto g_xx_0_0_0_x_y_x_yz = buffer_2000_pppd[22];

    auto g_xx_0_0_0_x_y_x_zz = buffer_2000_pppd[23];

    auto g_xx_0_0_0_x_y_y_xx = buffer_2000_pppd[24];

    auto g_xx_0_0_0_x_y_y_xy = buffer_2000_pppd[25];

    auto g_xx_0_0_0_x_y_y_xz = buffer_2000_pppd[26];

    auto g_xx_0_0_0_x_y_y_yy = buffer_2000_pppd[27];

    auto g_xx_0_0_0_x_y_y_yz = buffer_2000_pppd[28];

    auto g_xx_0_0_0_x_y_y_zz = buffer_2000_pppd[29];

    auto g_xx_0_0_0_x_y_z_xx = buffer_2000_pppd[30];

    auto g_xx_0_0_0_x_y_z_xy = buffer_2000_pppd[31];

    auto g_xx_0_0_0_x_y_z_xz = buffer_2000_pppd[32];

    auto g_xx_0_0_0_x_y_z_yy = buffer_2000_pppd[33];

    auto g_xx_0_0_0_x_y_z_yz = buffer_2000_pppd[34];

    auto g_xx_0_0_0_x_y_z_zz = buffer_2000_pppd[35];

    auto g_xx_0_0_0_x_z_x_xx = buffer_2000_pppd[36];

    auto g_xx_0_0_0_x_z_x_xy = buffer_2000_pppd[37];

    auto g_xx_0_0_0_x_z_x_xz = buffer_2000_pppd[38];

    auto g_xx_0_0_0_x_z_x_yy = buffer_2000_pppd[39];

    auto g_xx_0_0_0_x_z_x_yz = buffer_2000_pppd[40];

    auto g_xx_0_0_0_x_z_x_zz = buffer_2000_pppd[41];

    auto g_xx_0_0_0_x_z_y_xx = buffer_2000_pppd[42];

    auto g_xx_0_0_0_x_z_y_xy = buffer_2000_pppd[43];

    auto g_xx_0_0_0_x_z_y_xz = buffer_2000_pppd[44];

    auto g_xx_0_0_0_x_z_y_yy = buffer_2000_pppd[45];

    auto g_xx_0_0_0_x_z_y_yz = buffer_2000_pppd[46];

    auto g_xx_0_0_0_x_z_y_zz = buffer_2000_pppd[47];

    auto g_xx_0_0_0_x_z_z_xx = buffer_2000_pppd[48];

    auto g_xx_0_0_0_x_z_z_xy = buffer_2000_pppd[49];

    auto g_xx_0_0_0_x_z_z_xz = buffer_2000_pppd[50];

    auto g_xx_0_0_0_x_z_z_yy = buffer_2000_pppd[51];

    auto g_xx_0_0_0_x_z_z_yz = buffer_2000_pppd[52];

    auto g_xx_0_0_0_x_z_z_zz = buffer_2000_pppd[53];

    auto g_xx_0_0_0_y_x_x_xx = buffer_2000_pppd[54];

    auto g_xx_0_0_0_y_x_x_xy = buffer_2000_pppd[55];

    auto g_xx_0_0_0_y_x_x_xz = buffer_2000_pppd[56];

    auto g_xx_0_0_0_y_x_x_yy = buffer_2000_pppd[57];

    auto g_xx_0_0_0_y_x_x_yz = buffer_2000_pppd[58];

    auto g_xx_0_0_0_y_x_x_zz = buffer_2000_pppd[59];

    auto g_xx_0_0_0_y_x_y_xx = buffer_2000_pppd[60];

    auto g_xx_0_0_0_y_x_y_xy = buffer_2000_pppd[61];

    auto g_xx_0_0_0_y_x_y_xz = buffer_2000_pppd[62];

    auto g_xx_0_0_0_y_x_y_yy = buffer_2000_pppd[63];

    auto g_xx_0_0_0_y_x_y_yz = buffer_2000_pppd[64];

    auto g_xx_0_0_0_y_x_y_zz = buffer_2000_pppd[65];

    auto g_xx_0_0_0_y_x_z_xx = buffer_2000_pppd[66];

    auto g_xx_0_0_0_y_x_z_xy = buffer_2000_pppd[67];

    auto g_xx_0_0_0_y_x_z_xz = buffer_2000_pppd[68];

    auto g_xx_0_0_0_y_x_z_yy = buffer_2000_pppd[69];

    auto g_xx_0_0_0_y_x_z_yz = buffer_2000_pppd[70];

    auto g_xx_0_0_0_y_x_z_zz = buffer_2000_pppd[71];

    auto g_xx_0_0_0_y_y_x_xx = buffer_2000_pppd[72];

    auto g_xx_0_0_0_y_y_x_xy = buffer_2000_pppd[73];

    auto g_xx_0_0_0_y_y_x_xz = buffer_2000_pppd[74];

    auto g_xx_0_0_0_y_y_x_yy = buffer_2000_pppd[75];

    auto g_xx_0_0_0_y_y_x_yz = buffer_2000_pppd[76];

    auto g_xx_0_0_0_y_y_x_zz = buffer_2000_pppd[77];

    auto g_xx_0_0_0_y_y_y_xx = buffer_2000_pppd[78];

    auto g_xx_0_0_0_y_y_y_xy = buffer_2000_pppd[79];

    auto g_xx_0_0_0_y_y_y_xz = buffer_2000_pppd[80];

    auto g_xx_0_0_0_y_y_y_yy = buffer_2000_pppd[81];

    auto g_xx_0_0_0_y_y_y_yz = buffer_2000_pppd[82];

    auto g_xx_0_0_0_y_y_y_zz = buffer_2000_pppd[83];

    auto g_xx_0_0_0_y_y_z_xx = buffer_2000_pppd[84];

    auto g_xx_0_0_0_y_y_z_xy = buffer_2000_pppd[85];

    auto g_xx_0_0_0_y_y_z_xz = buffer_2000_pppd[86];

    auto g_xx_0_0_0_y_y_z_yy = buffer_2000_pppd[87];

    auto g_xx_0_0_0_y_y_z_yz = buffer_2000_pppd[88];

    auto g_xx_0_0_0_y_y_z_zz = buffer_2000_pppd[89];

    auto g_xx_0_0_0_y_z_x_xx = buffer_2000_pppd[90];

    auto g_xx_0_0_0_y_z_x_xy = buffer_2000_pppd[91];

    auto g_xx_0_0_0_y_z_x_xz = buffer_2000_pppd[92];

    auto g_xx_0_0_0_y_z_x_yy = buffer_2000_pppd[93];

    auto g_xx_0_0_0_y_z_x_yz = buffer_2000_pppd[94];

    auto g_xx_0_0_0_y_z_x_zz = buffer_2000_pppd[95];

    auto g_xx_0_0_0_y_z_y_xx = buffer_2000_pppd[96];

    auto g_xx_0_0_0_y_z_y_xy = buffer_2000_pppd[97];

    auto g_xx_0_0_0_y_z_y_xz = buffer_2000_pppd[98];

    auto g_xx_0_0_0_y_z_y_yy = buffer_2000_pppd[99];

    auto g_xx_0_0_0_y_z_y_yz = buffer_2000_pppd[100];

    auto g_xx_0_0_0_y_z_y_zz = buffer_2000_pppd[101];

    auto g_xx_0_0_0_y_z_z_xx = buffer_2000_pppd[102];

    auto g_xx_0_0_0_y_z_z_xy = buffer_2000_pppd[103];

    auto g_xx_0_0_0_y_z_z_xz = buffer_2000_pppd[104];

    auto g_xx_0_0_0_y_z_z_yy = buffer_2000_pppd[105];

    auto g_xx_0_0_0_y_z_z_yz = buffer_2000_pppd[106];

    auto g_xx_0_0_0_y_z_z_zz = buffer_2000_pppd[107];

    auto g_xx_0_0_0_z_x_x_xx = buffer_2000_pppd[108];

    auto g_xx_0_0_0_z_x_x_xy = buffer_2000_pppd[109];

    auto g_xx_0_0_0_z_x_x_xz = buffer_2000_pppd[110];

    auto g_xx_0_0_0_z_x_x_yy = buffer_2000_pppd[111];

    auto g_xx_0_0_0_z_x_x_yz = buffer_2000_pppd[112];

    auto g_xx_0_0_0_z_x_x_zz = buffer_2000_pppd[113];

    auto g_xx_0_0_0_z_x_y_xx = buffer_2000_pppd[114];

    auto g_xx_0_0_0_z_x_y_xy = buffer_2000_pppd[115];

    auto g_xx_0_0_0_z_x_y_xz = buffer_2000_pppd[116];

    auto g_xx_0_0_0_z_x_y_yy = buffer_2000_pppd[117];

    auto g_xx_0_0_0_z_x_y_yz = buffer_2000_pppd[118];

    auto g_xx_0_0_0_z_x_y_zz = buffer_2000_pppd[119];

    auto g_xx_0_0_0_z_x_z_xx = buffer_2000_pppd[120];

    auto g_xx_0_0_0_z_x_z_xy = buffer_2000_pppd[121];

    auto g_xx_0_0_0_z_x_z_xz = buffer_2000_pppd[122];

    auto g_xx_0_0_0_z_x_z_yy = buffer_2000_pppd[123];

    auto g_xx_0_0_0_z_x_z_yz = buffer_2000_pppd[124];

    auto g_xx_0_0_0_z_x_z_zz = buffer_2000_pppd[125];

    auto g_xx_0_0_0_z_y_x_xx = buffer_2000_pppd[126];

    auto g_xx_0_0_0_z_y_x_xy = buffer_2000_pppd[127];

    auto g_xx_0_0_0_z_y_x_xz = buffer_2000_pppd[128];

    auto g_xx_0_0_0_z_y_x_yy = buffer_2000_pppd[129];

    auto g_xx_0_0_0_z_y_x_yz = buffer_2000_pppd[130];

    auto g_xx_0_0_0_z_y_x_zz = buffer_2000_pppd[131];

    auto g_xx_0_0_0_z_y_y_xx = buffer_2000_pppd[132];

    auto g_xx_0_0_0_z_y_y_xy = buffer_2000_pppd[133];

    auto g_xx_0_0_0_z_y_y_xz = buffer_2000_pppd[134];

    auto g_xx_0_0_0_z_y_y_yy = buffer_2000_pppd[135];

    auto g_xx_0_0_0_z_y_y_yz = buffer_2000_pppd[136];

    auto g_xx_0_0_0_z_y_y_zz = buffer_2000_pppd[137];

    auto g_xx_0_0_0_z_y_z_xx = buffer_2000_pppd[138];

    auto g_xx_0_0_0_z_y_z_xy = buffer_2000_pppd[139];

    auto g_xx_0_0_0_z_y_z_xz = buffer_2000_pppd[140];

    auto g_xx_0_0_0_z_y_z_yy = buffer_2000_pppd[141];

    auto g_xx_0_0_0_z_y_z_yz = buffer_2000_pppd[142];

    auto g_xx_0_0_0_z_y_z_zz = buffer_2000_pppd[143];

    auto g_xx_0_0_0_z_z_x_xx = buffer_2000_pppd[144];

    auto g_xx_0_0_0_z_z_x_xy = buffer_2000_pppd[145];

    auto g_xx_0_0_0_z_z_x_xz = buffer_2000_pppd[146];

    auto g_xx_0_0_0_z_z_x_yy = buffer_2000_pppd[147];

    auto g_xx_0_0_0_z_z_x_yz = buffer_2000_pppd[148];

    auto g_xx_0_0_0_z_z_x_zz = buffer_2000_pppd[149];

    auto g_xx_0_0_0_z_z_y_xx = buffer_2000_pppd[150];

    auto g_xx_0_0_0_z_z_y_xy = buffer_2000_pppd[151];

    auto g_xx_0_0_0_z_z_y_xz = buffer_2000_pppd[152];

    auto g_xx_0_0_0_z_z_y_yy = buffer_2000_pppd[153];

    auto g_xx_0_0_0_z_z_y_yz = buffer_2000_pppd[154];

    auto g_xx_0_0_0_z_z_y_zz = buffer_2000_pppd[155];

    auto g_xx_0_0_0_z_z_z_xx = buffer_2000_pppd[156];

    auto g_xx_0_0_0_z_z_z_xy = buffer_2000_pppd[157];

    auto g_xx_0_0_0_z_z_z_xz = buffer_2000_pppd[158];

    auto g_xx_0_0_0_z_z_z_yy = buffer_2000_pppd[159];

    auto g_xx_0_0_0_z_z_z_yz = buffer_2000_pppd[160];

    auto g_xx_0_0_0_z_z_z_zz = buffer_2000_pppd[161];

    auto g_xy_0_0_0_x_x_x_xx = buffer_2000_pppd[162];

    auto g_xy_0_0_0_x_x_x_xy = buffer_2000_pppd[163];

    auto g_xy_0_0_0_x_x_x_xz = buffer_2000_pppd[164];

    auto g_xy_0_0_0_x_x_x_yy = buffer_2000_pppd[165];

    auto g_xy_0_0_0_x_x_x_yz = buffer_2000_pppd[166];

    auto g_xy_0_0_0_x_x_x_zz = buffer_2000_pppd[167];

    auto g_xy_0_0_0_x_x_y_xx = buffer_2000_pppd[168];

    auto g_xy_0_0_0_x_x_y_xy = buffer_2000_pppd[169];

    auto g_xy_0_0_0_x_x_y_xz = buffer_2000_pppd[170];

    auto g_xy_0_0_0_x_x_y_yy = buffer_2000_pppd[171];

    auto g_xy_0_0_0_x_x_y_yz = buffer_2000_pppd[172];

    auto g_xy_0_0_0_x_x_y_zz = buffer_2000_pppd[173];

    auto g_xy_0_0_0_x_x_z_xx = buffer_2000_pppd[174];

    auto g_xy_0_0_0_x_x_z_xy = buffer_2000_pppd[175];

    auto g_xy_0_0_0_x_x_z_xz = buffer_2000_pppd[176];

    auto g_xy_0_0_0_x_x_z_yy = buffer_2000_pppd[177];

    auto g_xy_0_0_0_x_x_z_yz = buffer_2000_pppd[178];

    auto g_xy_0_0_0_x_x_z_zz = buffer_2000_pppd[179];

    auto g_xy_0_0_0_x_y_x_xx = buffer_2000_pppd[180];

    auto g_xy_0_0_0_x_y_x_xy = buffer_2000_pppd[181];

    auto g_xy_0_0_0_x_y_x_xz = buffer_2000_pppd[182];

    auto g_xy_0_0_0_x_y_x_yy = buffer_2000_pppd[183];

    auto g_xy_0_0_0_x_y_x_yz = buffer_2000_pppd[184];

    auto g_xy_0_0_0_x_y_x_zz = buffer_2000_pppd[185];

    auto g_xy_0_0_0_x_y_y_xx = buffer_2000_pppd[186];

    auto g_xy_0_0_0_x_y_y_xy = buffer_2000_pppd[187];

    auto g_xy_0_0_0_x_y_y_xz = buffer_2000_pppd[188];

    auto g_xy_0_0_0_x_y_y_yy = buffer_2000_pppd[189];

    auto g_xy_0_0_0_x_y_y_yz = buffer_2000_pppd[190];

    auto g_xy_0_0_0_x_y_y_zz = buffer_2000_pppd[191];

    auto g_xy_0_0_0_x_y_z_xx = buffer_2000_pppd[192];

    auto g_xy_0_0_0_x_y_z_xy = buffer_2000_pppd[193];

    auto g_xy_0_0_0_x_y_z_xz = buffer_2000_pppd[194];

    auto g_xy_0_0_0_x_y_z_yy = buffer_2000_pppd[195];

    auto g_xy_0_0_0_x_y_z_yz = buffer_2000_pppd[196];

    auto g_xy_0_0_0_x_y_z_zz = buffer_2000_pppd[197];

    auto g_xy_0_0_0_x_z_x_xx = buffer_2000_pppd[198];

    auto g_xy_0_0_0_x_z_x_xy = buffer_2000_pppd[199];

    auto g_xy_0_0_0_x_z_x_xz = buffer_2000_pppd[200];

    auto g_xy_0_0_0_x_z_x_yy = buffer_2000_pppd[201];

    auto g_xy_0_0_0_x_z_x_yz = buffer_2000_pppd[202];

    auto g_xy_0_0_0_x_z_x_zz = buffer_2000_pppd[203];

    auto g_xy_0_0_0_x_z_y_xx = buffer_2000_pppd[204];

    auto g_xy_0_0_0_x_z_y_xy = buffer_2000_pppd[205];

    auto g_xy_0_0_0_x_z_y_xz = buffer_2000_pppd[206];

    auto g_xy_0_0_0_x_z_y_yy = buffer_2000_pppd[207];

    auto g_xy_0_0_0_x_z_y_yz = buffer_2000_pppd[208];

    auto g_xy_0_0_0_x_z_y_zz = buffer_2000_pppd[209];

    auto g_xy_0_0_0_x_z_z_xx = buffer_2000_pppd[210];

    auto g_xy_0_0_0_x_z_z_xy = buffer_2000_pppd[211];

    auto g_xy_0_0_0_x_z_z_xz = buffer_2000_pppd[212];

    auto g_xy_0_0_0_x_z_z_yy = buffer_2000_pppd[213];

    auto g_xy_0_0_0_x_z_z_yz = buffer_2000_pppd[214];

    auto g_xy_0_0_0_x_z_z_zz = buffer_2000_pppd[215];

    auto g_xy_0_0_0_y_x_x_xx = buffer_2000_pppd[216];

    auto g_xy_0_0_0_y_x_x_xy = buffer_2000_pppd[217];

    auto g_xy_0_0_0_y_x_x_xz = buffer_2000_pppd[218];

    auto g_xy_0_0_0_y_x_x_yy = buffer_2000_pppd[219];

    auto g_xy_0_0_0_y_x_x_yz = buffer_2000_pppd[220];

    auto g_xy_0_0_0_y_x_x_zz = buffer_2000_pppd[221];

    auto g_xy_0_0_0_y_x_y_xx = buffer_2000_pppd[222];

    auto g_xy_0_0_0_y_x_y_xy = buffer_2000_pppd[223];

    auto g_xy_0_0_0_y_x_y_xz = buffer_2000_pppd[224];

    auto g_xy_0_0_0_y_x_y_yy = buffer_2000_pppd[225];

    auto g_xy_0_0_0_y_x_y_yz = buffer_2000_pppd[226];

    auto g_xy_0_0_0_y_x_y_zz = buffer_2000_pppd[227];

    auto g_xy_0_0_0_y_x_z_xx = buffer_2000_pppd[228];

    auto g_xy_0_0_0_y_x_z_xy = buffer_2000_pppd[229];

    auto g_xy_0_0_0_y_x_z_xz = buffer_2000_pppd[230];

    auto g_xy_0_0_0_y_x_z_yy = buffer_2000_pppd[231];

    auto g_xy_0_0_0_y_x_z_yz = buffer_2000_pppd[232];

    auto g_xy_0_0_0_y_x_z_zz = buffer_2000_pppd[233];

    auto g_xy_0_0_0_y_y_x_xx = buffer_2000_pppd[234];

    auto g_xy_0_0_0_y_y_x_xy = buffer_2000_pppd[235];

    auto g_xy_0_0_0_y_y_x_xz = buffer_2000_pppd[236];

    auto g_xy_0_0_0_y_y_x_yy = buffer_2000_pppd[237];

    auto g_xy_0_0_0_y_y_x_yz = buffer_2000_pppd[238];

    auto g_xy_0_0_0_y_y_x_zz = buffer_2000_pppd[239];

    auto g_xy_0_0_0_y_y_y_xx = buffer_2000_pppd[240];

    auto g_xy_0_0_0_y_y_y_xy = buffer_2000_pppd[241];

    auto g_xy_0_0_0_y_y_y_xz = buffer_2000_pppd[242];

    auto g_xy_0_0_0_y_y_y_yy = buffer_2000_pppd[243];

    auto g_xy_0_0_0_y_y_y_yz = buffer_2000_pppd[244];

    auto g_xy_0_0_0_y_y_y_zz = buffer_2000_pppd[245];

    auto g_xy_0_0_0_y_y_z_xx = buffer_2000_pppd[246];

    auto g_xy_0_0_0_y_y_z_xy = buffer_2000_pppd[247];

    auto g_xy_0_0_0_y_y_z_xz = buffer_2000_pppd[248];

    auto g_xy_0_0_0_y_y_z_yy = buffer_2000_pppd[249];

    auto g_xy_0_0_0_y_y_z_yz = buffer_2000_pppd[250];

    auto g_xy_0_0_0_y_y_z_zz = buffer_2000_pppd[251];

    auto g_xy_0_0_0_y_z_x_xx = buffer_2000_pppd[252];

    auto g_xy_0_0_0_y_z_x_xy = buffer_2000_pppd[253];

    auto g_xy_0_0_0_y_z_x_xz = buffer_2000_pppd[254];

    auto g_xy_0_0_0_y_z_x_yy = buffer_2000_pppd[255];

    auto g_xy_0_0_0_y_z_x_yz = buffer_2000_pppd[256];

    auto g_xy_0_0_0_y_z_x_zz = buffer_2000_pppd[257];

    auto g_xy_0_0_0_y_z_y_xx = buffer_2000_pppd[258];

    auto g_xy_0_0_0_y_z_y_xy = buffer_2000_pppd[259];

    auto g_xy_0_0_0_y_z_y_xz = buffer_2000_pppd[260];

    auto g_xy_0_0_0_y_z_y_yy = buffer_2000_pppd[261];

    auto g_xy_0_0_0_y_z_y_yz = buffer_2000_pppd[262];

    auto g_xy_0_0_0_y_z_y_zz = buffer_2000_pppd[263];

    auto g_xy_0_0_0_y_z_z_xx = buffer_2000_pppd[264];

    auto g_xy_0_0_0_y_z_z_xy = buffer_2000_pppd[265];

    auto g_xy_0_0_0_y_z_z_xz = buffer_2000_pppd[266];

    auto g_xy_0_0_0_y_z_z_yy = buffer_2000_pppd[267];

    auto g_xy_0_0_0_y_z_z_yz = buffer_2000_pppd[268];

    auto g_xy_0_0_0_y_z_z_zz = buffer_2000_pppd[269];

    auto g_xy_0_0_0_z_x_x_xx = buffer_2000_pppd[270];

    auto g_xy_0_0_0_z_x_x_xy = buffer_2000_pppd[271];

    auto g_xy_0_0_0_z_x_x_xz = buffer_2000_pppd[272];

    auto g_xy_0_0_0_z_x_x_yy = buffer_2000_pppd[273];

    auto g_xy_0_0_0_z_x_x_yz = buffer_2000_pppd[274];

    auto g_xy_0_0_0_z_x_x_zz = buffer_2000_pppd[275];

    auto g_xy_0_0_0_z_x_y_xx = buffer_2000_pppd[276];

    auto g_xy_0_0_0_z_x_y_xy = buffer_2000_pppd[277];

    auto g_xy_0_0_0_z_x_y_xz = buffer_2000_pppd[278];

    auto g_xy_0_0_0_z_x_y_yy = buffer_2000_pppd[279];

    auto g_xy_0_0_0_z_x_y_yz = buffer_2000_pppd[280];

    auto g_xy_0_0_0_z_x_y_zz = buffer_2000_pppd[281];

    auto g_xy_0_0_0_z_x_z_xx = buffer_2000_pppd[282];

    auto g_xy_0_0_0_z_x_z_xy = buffer_2000_pppd[283];

    auto g_xy_0_0_0_z_x_z_xz = buffer_2000_pppd[284];

    auto g_xy_0_0_0_z_x_z_yy = buffer_2000_pppd[285];

    auto g_xy_0_0_0_z_x_z_yz = buffer_2000_pppd[286];

    auto g_xy_0_0_0_z_x_z_zz = buffer_2000_pppd[287];

    auto g_xy_0_0_0_z_y_x_xx = buffer_2000_pppd[288];

    auto g_xy_0_0_0_z_y_x_xy = buffer_2000_pppd[289];

    auto g_xy_0_0_0_z_y_x_xz = buffer_2000_pppd[290];

    auto g_xy_0_0_0_z_y_x_yy = buffer_2000_pppd[291];

    auto g_xy_0_0_0_z_y_x_yz = buffer_2000_pppd[292];

    auto g_xy_0_0_0_z_y_x_zz = buffer_2000_pppd[293];

    auto g_xy_0_0_0_z_y_y_xx = buffer_2000_pppd[294];

    auto g_xy_0_0_0_z_y_y_xy = buffer_2000_pppd[295];

    auto g_xy_0_0_0_z_y_y_xz = buffer_2000_pppd[296];

    auto g_xy_0_0_0_z_y_y_yy = buffer_2000_pppd[297];

    auto g_xy_0_0_0_z_y_y_yz = buffer_2000_pppd[298];

    auto g_xy_0_0_0_z_y_y_zz = buffer_2000_pppd[299];

    auto g_xy_0_0_0_z_y_z_xx = buffer_2000_pppd[300];

    auto g_xy_0_0_0_z_y_z_xy = buffer_2000_pppd[301];

    auto g_xy_0_0_0_z_y_z_xz = buffer_2000_pppd[302];

    auto g_xy_0_0_0_z_y_z_yy = buffer_2000_pppd[303];

    auto g_xy_0_0_0_z_y_z_yz = buffer_2000_pppd[304];

    auto g_xy_0_0_0_z_y_z_zz = buffer_2000_pppd[305];

    auto g_xy_0_0_0_z_z_x_xx = buffer_2000_pppd[306];

    auto g_xy_0_0_0_z_z_x_xy = buffer_2000_pppd[307];

    auto g_xy_0_0_0_z_z_x_xz = buffer_2000_pppd[308];

    auto g_xy_0_0_0_z_z_x_yy = buffer_2000_pppd[309];

    auto g_xy_0_0_0_z_z_x_yz = buffer_2000_pppd[310];

    auto g_xy_0_0_0_z_z_x_zz = buffer_2000_pppd[311];

    auto g_xy_0_0_0_z_z_y_xx = buffer_2000_pppd[312];

    auto g_xy_0_0_0_z_z_y_xy = buffer_2000_pppd[313];

    auto g_xy_0_0_0_z_z_y_xz = buffer_2000_pppd[314];

    auto g_xy_0_0_0_z_z_y_yy = buffer_2000_pppd[315];

    auto g_xy_0_0_0_z_z_y_yz = buffer_2000_pppd[316];

    auto g_xy_0_0_0_z_z_y_zz = buffer_2000_pppd[317];

    auto g_xy_0_0_0_z_z_z_xx = buffer_2000_pppd[318];

    auto g_xy_0_0_0_z_z_z_xy = buffer_2000_pppd[319];

    auto g_xy_0_0_0_z_z_z_xz = buffer_2000_pppd[320];

    auto g_xy_0_0_0_z_z_z_yy = buffer_2000_pppd[321];

    auto g_xy_0_0_0_z_z_z_yz = buffer_2000_pppd[322];

    auto g_xy_0_0_0_z_z_z_zz = buffer_2000_pppd[323];

    auto g_xz_0_0_0_x_x_x_xx = buffer_2000_pppd[324];

    auto g_xz_0_0_0_x_x_x_xy = buffer_2000_pppd[325];

    auto g_xz_0_0_0_x_x_x_xz = buffer_2000_pppd[326];

    auto g_xz_0_0_0_x_x_x_yy = buffer_2000_pppd[327];

    auto g_xz_0_0_0_x_x_x_yz = buffer_2000_pppd[328];

    auto g_xz_0_0_0_x_x_x_zz = buffer_2000_pppd[329];

    auto g_xz_0_0_0_x_x_y_xx = buffer_2000_pppd[330];

    auto g_xz_0_0_0_x_x_y_xy = buffer_2000_pppd[331];

    auto g_xz_0_0_0_x_x_y_xz = buffer_2000_pppd[332];

    auto g_xz_0_0_0_x_x_y_yy = buffer_2000_pppd[333];

    auto g_xz_0_0_0_x_x_y_yz = buffer_2000_pppd[334];

    auto g_xz_0_0_0_x_x_y_zz = buffer_2000_pppd[335];

    auto g_xz_0_0_0_x_x_z_xx = buffer_2000_pppd[336];

    auto g_xz_0_0_0_x_x_z_xy = buffer_2000_pppd[337];

    auto g_xz_0_0_0_x_x_z_xz = buffer_2000_pppd[338];

    auto g_xz_0_0_0_x_x_z_yy = buffer_2000_pppd[339];

    auto g_xz_0_0_0_x_x_z_yz = buffer_2000_pppd[340];

    auto g_xz_0_0_0_x_x_z_zz = buffer_2000_pppd[341];

    auto g_xz_0_0_0_x_y_x_xx = buffer_2000_pppd[342];

    auto g_xz_0_0_0_x_y_x_xy = buffer_2000_pppd[343];

    auto g_xz_0_0_0_x_y_x_xz = buffer_2000_pppd[344];

    auto g_xz_0_0_0_x_y_x_yy = buffer_2000_pppd[345];

    auto g_xz_0_0_0_x_y_x_yz = buffer_2000_pppd[346];

    auto g_xz_0_0_0_x_y_x_zz = buffer_2000_pppd[347];

    auto g_xz_0_0_0_x_y_y_xx = buffer_2000_pppd[348];

    auto g_xz_0_0_0_x_y_y_xy = buffer_2000_pppd[349];

    auto g_xz_0_0_0_x_y_y_xz = buffer_2000_pppd[350];

    auto g_xz_0_0_0_x_y_y_yy = buffer_2000_pppd[351];

    auto g_xz_0_0_0_x_y_y_yz = buffer_2000_pppd[352];

    auto g_xz_0_0_0_x_y_y_zz = buffer_2000_pppd[353];

    auto g_xz_0_0_0_x_y_z_xx = buffer_2000_pppd[354];

    auto g_xz_0_0_0_x_y_z_xy = buffer_2000_pppd[355];

    auto g_xz_0_0_0_x_y_z_xz = buffer_2000_pppd[356];

    auto g_xz_0_0_0_x_y_z_yy = buffer_2000_pppd[357];

    auto g_xz_0_0_0_x_y_z_yz = buffer_2000_pppd[358];

    auto g_xz_0_0_0_x_y_z_zz = buffer_2000_pppd[359];

    auto g_xz_0_0_0_x_z_x_xx = buffer_2000_pppd[360];

    auto g_xz_0_0_0_x_z_x_xy = buffer_2000_pppd[361];

    auto g_xz_0_0_0_x_z_x_xz = buffer_2000_pppd[362];

    auto g_xz_0_0_0_x_z_x_yy = buffer_2000_pppd[363];

    auto g_xz_0_0_0_x_z_x_yz = buffer_2000_pppd[364];

    auto g_xz_0_0_0_x_z_x_zz = buffer_2000_pppd[365];

    auto g_xz_0_0_0_x_z_y_xx = buffer_2000_pppd[366];

    auto g_xz_0_0_0_x_z_y_xy = buffer_2000_pppd[367];

    auto g_xz_0_0_0_x_z_y_xz = buffer_2000_pppd[368];

    auto g_xz_0_0_0_x_z_y_yy = buffer_2000_pppd[369];

    auto g_xz_0_0_0_x_z_y_yz = buffer_2000_pppd[370];

    auto g_xz_0_0_0_x_z_y_zz = buffer_2000_pppd[371];

    auto g_xz_0_0_0_x_z_z_xx = buffer_2000_pppd[372];

    auto g_xz_0_0_0_x_z_z_xy = buffer_2000_pppd[373];

    auto g_xz_0_0_0_x_z_z_xz = buffer_2000_pppd[374];

    auto g_xz_0_0_0_x_z_z_yy = buffer_2000_pppd[375];

    auto g_xz_0_0_0_x_z_z_yz = buffer_2000_pppd[376];

    auto g_xz_0_0_0_x_z_z_zz = buffer_2000_pppd[377];

    auto g_xz_0_0_0_y_x_x_xx = buffer_2000_pppd[378];

    auto g_xz_0_0_0_y_x_x_xy = buffer_2000_pppd[379];

    auto g_xz_0_0_0_y_x_x_xz = buffer_2000_pppd[380];

    auto g_xz_0_0_0_y_x_x_yy = buffer_2000_pppd[381];

    auto g_xz_0_0_0_y_x_x_yz = buffer_2000_pppd[382];

    auto g_xz_0_0_0_y_x_x_zz = buffer_2000_pppd[383];

    auto g_xz_0_0_0_y_x_y_xx = buffer_2000_pppd[384];

    auto g_xz_0_0_0_y_x_y_xy = buffer_2000_pppd[385];

    auto g_xz_0_0_0_y_x_y_xz = buffer_2000_pppd[386];

    auto g_xz_0_0_0_y_x_y_yy = buffer_2000_pppd[387];

    auto g_xz_0_0_0_y_x_y_yz = buffer_2000_pppd[388];

    auto g_xz_0_0_0_y_x_y_zz = buffer_2000_pppd[389];

    auto g_xz_0_0_0_y_x_z_xx = buffer_2000_pppd[390];

    auto g_xz_0_0_0_y_x_z_xy = buffer_2000_pppd[391];

    auto g_xz_0_0_0_y_x_z_xz = buffer_2000_pppd[392];

    auto g_xz_0_0_0_y_x_z_yy = buffer_2000_pppd[393];

    auto g_xz_0_0_0_y_x_z_yz = buffer_2000_pppd[394];

    auto g_xz_0_0_0_y_x_z_zz = buffer_2000_pppd[395];

    auto g_xz_0_0_0_y_y_x_xx = buffer_2000_pppd[396];

    auto g_xz_0_0_0_y_y_x_xy = buffer_2000_pppd[397];

    auto g_xz_0_0_0_y_y_x_xz = buffer_2000_pppd[398];

    auto g_xz_0_0_0_y_y_x_yy = buffer_2000_pppd[399];

    auto g_xz_0_0_0_y_y_x_yz = buffer_2000_pppd[400];

    auto g_xz_0_0_0_y_y_x_zz = buffer_2000_pppd[401];

    auto g_xz_0_0_0_y_y_y_xx = buffer_2000_pppd[402];

    auto g_xz_0_0_0_y_y_y_xy = buffer_2000_pppd[403];

    auto g_xz_0_0_0_y_y_y_xz = buffer_2000_pppd[404];

    auto g_xz_0_0_0_y_y_y_yy = buffer_2000_pppd[405];

    auto g_xz_0_0_0_y_y_y_yz = buffer_2000_pppd[406];

    auto g_xz_0_0_0_y_y_y_zz = buffer_2000_pppd[407];

    auto g_xz_0_0_0_y_y_z_xx = buffer_2000_pppd[408];

    auto g_xz_0_0_0_y_y_z_xy = buffer_2000_pppd[409];

    auto g_xz_0_0_0_y_y_z_xz = buffer_2000_pppd[410];

    auto g_xz_0_0_0_y_y_z_yy = buffer_2000_pppd[411];

    auto g_xz_0_0_0_y_y_z_yz = buffer_2000_pppd[412];

    auto g_xz_0_0_0_y_y_z_zz = buffer_2000_pppd[413];

    auto g_xz_0_0_0_y_z_x_xx = buffer_2000_pppd[414];

    auto g_xz_0_0_0_y_z_x_xy = buffer_2000_pppd[415];

    auto g_xz_0_0_0_y_z_x_xz = buffer_2000_pppd[416];

    auto g_xz_0_0_0_y_z_x_yy = buffer_2000_pppd[417];

    auto g_xz_0_0_0_y_z_x_yz = buffer_2000_pppd[418];

    auto g_xz_0_0_0_y_z_x_zz = buffer_2000_pppd[419];

    auto g_xz_0_0_0_y_z_y_xx = buffer_2000_pppd[420];

    auto g_xz_0_0_0_y_z_y_xy = buffer_2000_pppd[421];

    auto g_xz_0_0_0_y_z_y_xz = buffer_2000_pppd[422];

    auto g_xz_0_0_0_y_z_y_yy = buffer_2000_pppd[423];

    auto g_xz_0_0_0_y_z_y_yz = buffer_2000_pppd[424];

    auto g_xz_0_0_0_y_z_y_zz = buffer_2000_pppd[425];

    auto g_xz_0_0_0_y_z_z_xx = buffer_2000_pppd[426];

    auto g_xz_0_0_0_y_z_z_xy = buffer_2000_pppd[427];

    auto g_xz_0_0_0_y_z_z_xz = buffer_2000_pppd[428];

    auto g_xz_0_0_0_y_z_z_yy = buffer_2000_pppd[429];

    auto g_xz_0_0_0_y_z_z_yz = buffer_2000_pppd[430];

    auto g_xz_0_0_0_y_z_z_zz = buffer_2000_pppd[431];

    auto g_xz_0_0_0_z_x_x_xx = buffer_2000_pppd[432];

    auto g_xz_0_0_0_z_x_x_xy = buffer_2000_pppd[433];

    auto g_xz_0_0_0_z_x_x_xz = buffer_2000_pppd[434];

    auto g_xz_0_0_0_z_x_x_yy = buffer_2000_pppd[435];

    auto g_xz_0_0_0_z_x_x_yz = buffer_2000_pppd[436];

    auto g_xz_0_0_0_z_x_x_zz = buffer_2000_pppd[437];

    auto g_xz_0_0_0_z_x_y_xx = buffer_2000_pppd[438];

    auto g_xz_0_0_0_z_x_y_xy = buffer_2000_pppd[439];

    auto g_xz_0_0_0_z_x_y_xz = buffer_2000_pppd[440];

    auto g_xz_0_0_0_z_x_y_yy = buffer_2000_pppd[441];

    auto g_xz_0_0_0_z_x_y_yz = buffer_2000_pppd[442];

    auto g_xz_0_0_0_z_x_y_zz = buffer_2000_pppd[443];

    auto g_xz_0_0_0_z_x_z_xx = buffer_2000_pppd[444];

    auto g_xz_0_0_0_z_x_z_xy = buffer_2000_pppd[445];

    auto g_xz_0_0_0_z_x_z_xz = buffer_2000_pppd[446];

    auto g_xz_0_0_0_z_x_z_yy = buffer_2000_pppd[447];

    auto g_xz_0_0_0_z_x_z_yz = buffer_2000_pppd[448];

    auto g_xz_0_0_0_z_x_z_zz = buffer_2000_pppd[449];

    auto g_xz_0_0_0_z_y_x_xx = buffer_2000_pppd[450];

    auto g_xz_0_0_0_z_y_x_xy = buffer_2000_pppd[451];

    auto g_xz_0_0_0_z_y_x_xz = buffer_2000_pppd[452];

    auto g_xz_0_0_0_z_y_x_yy = buffer_2000_pppd[453];

    auto g_xz_0_0_0_z_y_x_yz = buffer_2000_pppd[454];

    auto g_xz_0_0_0_z_y_x_zz = buffer_2000_pppd[455];

    auto g_xz_0_0_0_z_y_y_xx = buffer_2000_pppd[456];

    auto g_xz_0_0_0_z_y_y_xy = buffer_2000_pppd[457];

    auto g_xz_0_0_0_z_y_y_xz = buffer_2000_pppd[458];

    auto g_xz_0_0_0_z_y_y_yy = buffer_2000_pppd[459];

    auto g_xz_0_0_0_z_y_y_yz = buffer_2000_pppd[460];

    auto g_xz_0_0_0_z_y_y_zz = buffer_2000_pppd[461];

    auto g_xz_0_0_0_z_y_z_xx = buffer_2000_pppd[462];

    auto g_xz_0_0_0_z_y_z_xy = buffer_2000_pppd[463];

    auto g_xz_0_0_0_z_y_z_xz = buffer_2000_pppd[464];

    auto g_xz_0_0_0_z_y_z_yy = buffer_2000_pppd[465];

    auto g_xz_0_0_0_z_y_z_yz = buffer_2000_pppd[466];

    auto g_xz_0_0_0_z_y_z_zz = buffer_2000_pppd[467];

    auto g_xz_0_0_0_z_z_x_xx = buffer_2000_pppd[468];

    auto g_xz_0_0_0_z_z_x_xy = buffer_2000_pppd[469];

    auto g_xz_0_0_0_z_z_x_xz = buffer_2000_pppd[470];

    auto g_xz_0_0_0_z_z_x_yy = buffer_2000_pppd[471];

    auto g_xz_0_0_0_z_z_x_yz = buffer_2000_pppd[472];

    auto g_xz_0_0_0_z_z_x_zz = buffer_2000_pppd[473];

    auto g_xz_0_0_0_z_z_y_xx = buffer_2000_pppd[474];

    auto g_xz_0_0_0_z_z_y_xy = buffer_2000_pppd[475];

    auto g_xz_0_0_0_z_z_y_xz = buffer_2000_pppd[476];

    auto g_xz_0_0_0_z_z_y_yy = buffer_2000_pppd[477];

    auto g_xz_0_0_0_z_z_y_yz = buffer_2000_pppd[478];

    auto g_xz_0_0_0_z_z_y_zz = buffer_2000_pppd[479];

    auto g_xz_0_0_0_z_z_z_xx = buffer_2000_pppd[480];

    auto g_xz_0_0_0_z_z_z_xy = buffer_2000_pppd[481];

    auto g_xz_0_0_0_z_z_z_xz = buffer_2000_pppd[482];

    auto g_xz_0_0_0_z_z_z_yy = buffer_2000_pppd[483];

    auto g_xz_0_0_0_z_z_z_yz = buffer_2000_pppd[484];

    auto g_xz_0_0_0_z_z_z_zz = buffer_2000_pppd[485];

    auto g_yy_0_0_0_x_x_x_xx = buffer_2000_pppd[486];

    auto g_yy_0_0_0_x_x_x_xy = buffer_2000_pppd[487];

    auto g_yy_0_0_0_x_x_x_xz = buffer_2000_pppd[488];

    auto g_yy_0_0_0_x_x_x_yy = buffer_2000_pppd[489];

    auto g_yy_0_0_0_x_x_x_yz = buffer_2000_pppd[490];

    auto g_yy_0_0_0_x_x_x_zz = buffer_2000_pppd[491];

    auto g_yy_0_0_0_x_x_y_xx = buffer_2000_pppd[492];

    auto g_yy_0_0_0_x_x_y_xy = buffer_2000_pppd[493];

    auto g_yy_0_0_0_x_x_y_xz = buffer_2000_pppd[494];

    auto g_yy_0_0_0_x_x_y_yy = buffer_2000_pppd[495];

    auto g_yy_0_0_0_x_x_y_yz = buffer_2000_pppd[496];

    auto g_yy_0_0_0_x_x_y_zz = buffer_2000_pppd[497];

    auto g_yy_0_0_0_x_x_z_xx = buffer_2000_pppd[498];

    auto g_yy_0_0_0_x_x_z_xy = buffer_2000_pppd[499];

    auto g_yy_0_0_0_x_x_z_xz = buffer_2000_pppd[500];

    auto g_yy_0_0_0_x_x_z_yy = buffer_2000_pppd[501];

    auto g_yy_0_0_0_x_x_z_yz = buffer_2000_pppd[502];

    auto g_yy_0_0_0_x_x_z_zz = buffer_2000_pppd[503];

    auto g_yy_0_0_0_x_y_x_xx = buffer_2000_pppd[504];

    auto g_yy_0_0_0_x_y_x_xy = buffer_2000_pppd[505];

    auto g_yy_0_0_0_x_y_x_xz = buffer_2000_pppd[506];

    auto g_yy_0_0_0_x_y_x_yy = buffer_2000_pppd[507];

    auto g_yy_0_0_0_x_y_x_yz = buffer_2000_pppd[508];

    auto g_yy_0_0_0_x_y_x_zz = buffer_2000_pppd[509];

    auto g_yy_0_0_0_x_y_y_xx = buffer_2000_pppd[510];

    auto g_yy_0_0_0_x_y_y_xy = buffer_2000_pppd[511];

    auto g_yy_0_0_0_x_y_y_xz = buffer_2000_pppd[512];

    auto g_yy_0_0_0_x_y_y_yy = buffer_2000_pppd[513];

    auto g_yy_0_0_0_x_y_y_yz = buffer_2000_pppd[514];

    auto g_yy_0_0_0_x_y_y_zz = buffer_2000_pppd[515];

    auto g_yy_0_0_0_x_y_z_xx = buffer_2000_pppd[516];

    auto g_yy_0_0_0_x_y_z_xy = buffer_2000_pppd[517];

    auto g_yy_0_0_0_x_y_z_xz = buffer_2000_pppd[518];

    auto g_yy_0_0_0_x_y_z_yy = buffer_2000_pppd[519];

    auto g_yy_0_0_0_x_y_z_yz = buffer_2000_pppd[520];

    auto g_yy_0_0_0_x_y_z_zz = buffer_2000_pppd[521];

    auto g_yy_0_0_0_x_z_x_xx = buffer_2000_pppd[522];

    auto g_yy_0_0_0_x_z_x_xy = buffer_2000_pppd[523];

    auto g_yy_0_0_0_x_z_x_xz = buffer_2000_pppd[524];

    auto g_yy_0_0_0_x_z_x_yy = buffer_2000_pppd[525];

    auto g_yy_0_0_0_x_z_x_yz = buffer_2000_pppd[526];

    auto g_yy_0_0_0_x_z_x_zz = buffer_2000_pppd[527];

    auto g_yy_0_0_0_x_z_y_xx = buffer_2000_pppd[528];

    auto g_yy_0_0_0_x_z_y_xy = buffer_2000_pppd[529];

    auto g_yy_0_0_0_x_z_y_xz = buffer_2000_pppd[530];

    auto g_yy_0_0_0_x_z_y_yy = buffer_2000_pppd[531];

    auto g_yy_0_0_0_x_z_y_yz = buffer_2000_pppd[532];

    auto g_yy_0_0_0_x_z_y_zz = buffer_2000_pppd[533];

    auto g_yy_0_0_0_x_z_z_xx = buffer_2000_pppd[534];

    auto g_yy_0_0_0_x_z_z_xy = buffer_2000_pppd[535];

    auto g_yy_0_0_0_x_z_z_xz = buffer_2000_pppd[536];

    auto g_yy_0_0_0_x_z_z_yy = buffer_2000_pppd[537];

    auto g_yy_0_0_0_x_z_z_yz = buffer_2000_pppd[538];

    auto g_yy_0_0_0_x_z_z_zz = buffer_2000_pppd[539];

    auto g_yy_0_0_0_y_x_x_xx = buffer_2000_pppd[540];

    auto g_yy_0_0_0_y_x_x_xy = buffer_2000_pppd[541];

    auto g_yy_0_0_0_y_x_x_xz = buffer_2000_pppd[542];

    auto g_yy_0_0_0_y_x_x_yy = buffer_2000_pppd[543];

    auto g_yy_0_0_0_y_x_x_yz = buffer_2000_pppd[544];

    auto g_yy_0_0_0_y_x_x_zz = buffer_2000_pppd[545];

    auto g_yy_0_0_0_y_x_y_xx = buffer_2000_pppd[546];

    auto g_yy_0_0_0_y_x_y_xy = buffer_2000_pppd[547];

    auto g_yy_0_0_0_y_x_y_xz = buffer_2000_pppd[548];

    auto g_yy_0_0_0_y_x_y_yy = buffer_2000_pppd[549];

    auto g_yy_0_0_0_y_x_y_yz = buffer_2000_pppd[550];

    auto g_yy_0_0_0_y_x_y_zz = buffer_2000_pppd[551];

    auto g_yy_0_0_0_y_x_z_xx = buffer_2000_pppd[552];

    auto g_yy_0_0_0_y_x_z_xy = buffer_2000_pppd[553];

    auto g_yy_0_0_0_y_x_z_xz = buffer_2000_pppd[554];

    auto g_yy_0_0_0_y_x_z_yy = buffer_2000_pppd[555];

    auto g_yy_0_0_0_y_x_z_yz = buffer_2000_pppd[556];

    auto g_yy_0_0_0_y_x_z_zz = buffer_2000_pppd[557];

    auto g_yy_0_0_0_y_y_x_xx = buffer_2000_pppd[558];

    auto g_yy_0_0_0_y_y_x_xy = buffer_2000_pppd[559];

    auto g_yy_0_0_0_y_y_x_xz = buffer_2000_pppd[560];

    auto g_yy_0_0_0_y_y_x_yy = buffer_2000_pppd[561];

    auto g_yy_0_0_0_y_y_x_yz = buffer_2000_pppd[562];

    auto g_yy_0_0_0_y_y_x_zz = buffer_2000_pppd[563];

    auto g_yy_0_0_0_y_y_y_xx = buffer_2000_pppd[564];

    auto g_yy_0_0_0_y_y_y_xy = buffer_2000_pppd[565];

    auto g_yy_0_0_0_y_y_y_xz = buffer_2000_pppd[566];

    auto g_yy_0_0_0_y_y_y_yy = buffer_2000_pppd[567];

    auto g_yy_0_0_0_y_y_y_yz = buffer_2000_pppd[568];

    auto g_yy_0_0_0_y_y_y_zz = buffer_2000_pppd[569];

    auto g_yy_0_0_0_y_y_z_xx = buffer_2000_pppd[570];

    auto g_yy_0_0_0_y_y_z_xy = buffer_2000_pppd[571];

    auto g_yy_0_0_0_y_y_z_xz = buffer_2000_pppd[572];

    auto g_yy_0_0_0_y_y_z_yy = buffer_2000_pppd[573];

    auto g_yy_0_0_0_y_y_z_yz = buffer_2000_pppd[574];

    auto g_yy_0_0_0_y_y_z_zz = buffer_2000_pppd[575];

    auto g_yy_0_0_0_y_z_x_xx = buffer_2000_pppd[576];

    auto g_yy_0_0_0_y_z_x_xy = buffer_2000_pppd[577];

    auto g_yy_0_0_0_y_z_x_xz = buffer_2000_pppd[578];

    auto g_yy_0_0_0_y_z_x_yy = buffer_2000_pppd[579];

    auto g_yy_0_0_0_y_z_x_yz = buffer_2000_pppd[580];

    auto g_yy_0_0_0_y_z_x_zz = buffer_2000_pppd[581];

    auto g_yy_0_0_0_y_z_y_xx = buffer_2000_pppd[582];

    auto g_yy_0_0_0_y_z_y_xy = buffer_2000_pppd[583];

    auto g_yy_0_0_0_y_z_y_xz = buffer_2000_pppd[584];

    auto g_yy_0_0_0_y_z_y_yy = buffer_2000_pppd[585];

    auto g_yy_0_0_0_y_z_y_yz = buffer_2000_pppd[586];

    auto g_yy_0_0_0_y_z_y_zz = buffer_2000_pppd[587];

    auto g_yy_0_0_0_y_z_z_xx = buffer_2000_pppd[588];

    auto g_yy_0_0_0_y_z_z_xy = buffer_2000_pppd[589];

    auto g_yy_0_0_0_y_z_z_xz = buffer_2000_pppd[590];

    auto g_yy_0_0_0_y_z_z_yy = buffer_2000_pppd[591];

    auto g_yy_0_0_0_y_z_z_yz = buffer_2000_pppd[592];

    auto g_yy_0_0_0_y_z_z_zz = buffer_2000_pppd[593];

    auto g_yy_0_0_0_z_x_x_xx = buffer_2000_pppd[594];

    auto g_yy_0_0_0_z_x_x_xy = buffer_2000_pppd[595];

    auto g_yy_0_0_0_z_x_x_xz = buffer_2000_pppd[596];

    auto g_yy_0_0_0_z_x_x_yy = buffer_2000_pppd[597];

    auto g_yy_0_0_0_z_x_x_yz = buffer_2000_pppd[598];

    auto g_yy_0_0_0_z_x_x_zz = buffer_2000_pppd[599];

    auto g_yy_0_0_0_z_x_y_xx = buffer_2000_pppd[600];

    auto g_yy_0_0_0_z_x_y_xy = buffer_2000_pppd[601];

    auto g_yy_0_0_0_z_x_y_xz = buffer_2000_pppd[602];

    auto g_yy_0_0_0_z_x_y_yy = buffer_2000_pppd[603];

    auto g_yy_0_0_0_z_x_y_yz = buffer_2000_pppd[604];

    auto g_yy_0_0_0_z_x_y_zz = buffer_2000_pppd[605];

    auto g_yy_0_0_0_z_x_z_xx = buffer_2000_pppd[606];

    auto g_yy_0_0_0_z_x_z_xy = buffer_2000_pppd[607];

    auto g_yy_0_0_0_z_x_z_xz = buffer_2000_pppd[608];

    auto g_yy_0_0_0_z_x_z_yy = buffer_2000_pppd[609];

    auto g_yy_0_0_0_z_x_z_yz = buffer_2000_pppd[610];

    auto g_yy_0_0_0_z_x_z_zz = buffer_2000_pppd[611];

    auto g_yy_0_0_0_z_y_x_xx = buffer_2000_pppd[612];

    auto g_yy_0_0_0_z_y_x_xy = buffer_2000_pppd[613];

    auto g_yy_0_0_0_z_y_x_xz = buffer_2000_pppd[614];

    auto g_yy_0_0_0_z_y_x_yy = buffer_2000_pppd[615];

    auto g_yy_0_0_0_z_y_x_yz = buffer_2000_pppd[616];

    auto g_yy_0_0_0_z_y_x_zz = buffer_2000_pppd[617];

    auto g_yy_0_0_0_z_y_y_xx = buffer_2000_pppd[618];

    auto g_yy_0_0_0_z_y_y_xy = buffer_2000_pppd[619];

    auto g_yy_0_0_0_z_y_y_xz = buffer_2000_pppd[620];

    auto g_yy_0_0_0_z_y_y_yy = buffer_2000_pppd[621];

    auto g_yy_0_0_0_z_y_y_yz = buffer_2000_pppd[622];

    auto g_yy_0_0_0_z_y_y_zz = buffer_2000_pppd[623];

    auto g_yy_0_0_0_z_y_z_xx = buffer_2000_pppd[624];

    auto g_yy_0_0_0_z_y_z_xy = buffer_2000_pppd[625];

    auto g_yy_0_0_0_z_y_z_xz = buffer_2000_pppd[626];

    auto g_yy_0_0_0_z_y_z_yy = buffer_2000_pppd[627];

    auto g_yy_0_0_0_z_y_z_yz = buffer_2000_pppd[628];

    auto g_yy_0_0_0_z_y_z_zz = buffer_2000_pppd[629];

    auto g_yy_0_0_0_z_z_x_xx = buffer_2000_pppd[630];

    auto g_yy_0_0_0_z_z_x_xy = buffer_2000_pppd[631];

    auto g_yy_0_0_0_z_z_x_xz = buffer_2000_pppd[632];

    auto g_yy_0_0_0_z_z_x_yy = buffer_2000_pppd[633];

    auto g_yy_0_0_0_z_z_x_yz = buffer_2000_pppd[634];

    auto g_yy_0_0_0_z_z_x_zz = buffer_2000_pppd[635];

    auto g_yy_0_0_0_z_z_y_xx = buffer_2000_pppd[636];

    auto g_yy_0_0_0_z_z_y_xy = buffer_2000_pppd[637];

    auto g_yy_0_0_0_z_z_y_xz = buffer_2000_pppd[638];

    auto g_yy_0_0_0_z_z_y_yy = buffer_2000_pppd[639];

    auto g_yy_0_0_0_z_z_y_yz = buffer_2000_pppd[640];

    auto g_yy_0_0_0_z_z_y_zz = buffer_2000_pppd[641];

    auto g_yy_0_0_0_z_z_z_xx = buffer_2000_pppd[642];

    auto g_yy_0_0_0_z_z_z_xy = buffer_2000_pppd[643];

    auto g_yy_0_0_0_z_z_z_xz = buffer_2000_pppd[644];

    auto g_yy_0_0_0_z_z_z_yy = buffer_2000_pppd[645];

    auto g_yy_0_0_0_z_z_z_yz = buffer_2000_pppd[646];

    auto g_yy_0_0_0_z_z_z_zz = buffer_2000_pppd[647];

    auto g_yz_0_0_0_x_x_x_xx = buffer_2000_pppd[648];

    auto g_yz_0_0_0_x_x_x_xy = buffer_2000_pppd[649];

    auto g_yz_0_0_0_x_x_x_xz = buffer_2000_pppd[650];

    auto g_yz_0_0_0_x_x_x_yy = buffer_2000_pppd[651];

    auto g_yz_0_0_0_x_x_x_yz = buffer_2000_pppd[652];

    auto g_yz_0_0_0_x_x_x_zz = buffer_2000_pppd[653];

    auto g_yz_0_0_0_x_x_y_xx = buffer_2000_pppd[654];

    auto g_yz_0_0_0_x_x_y_xy = buffer_2000_pppd[655];

    auto g_yz_0_0_0_x_x_y_xz = buffer_2000_pppd[656];

    auto g_yz_0_0_0_x_x_y_yy = buffer_2000_pppd[657];

    auto g_yz_0_0_0_x_x_y_yz = buffer_2000_pppd[658];

    auto g_yz_0_0_0_x_x_y_zz = buffer_2000_pppd[659];

    auto g_yz_0_0_0_x_x_z_xx = buffer_2000_pppd[660];

    auto g_yz_0_0_0_x_x_z_xy = buffer_2000_pppd[661];

    auto g_yz_0_0_0_x_x_z_xz = buffer_2000_pppd[662];

    auto g_yz_0_0_0_x_x_z_yy = buffer_2000_pppd[663];

    auto g_yz_0_0_0_x_x_z_yz = buffer_2000_pppd[664];

    auto g_yz_0_0_0_x_x_z_zz = buffer_2000_pppd[665];

    auto g_yz_0_0_0_x_y_x_xx = buffer_2000_pppd[666];

    auto g_yz_0_0_0_x_y_x_xy = buffer_2000_pppd[667];

    auto g_yz_0_0_0_x_y_x_xz = buffer_2000_pppd[668];

    auto g_yz_0_0_0_x_y_x_yy = buffer_2000_pppd[669];

    auto g_yz_0_0_0_x_y_x_yz = buffer_2000_pppd[670];

    auto g_yz_0_0_0_x_y_x_zz = buffer_2000_pppd[671];

    auto g_yz_0_0_0_x_y_y_xx = buffer_2000_pppd[672];

    auto g_yz_0_0_0_x_y_y_xy = buffer_2000_pppd[673];

    auto g_yz_0_0_0_x_y_y_xz = buffer_2000_pppd[674];

    auto g_yz_0_0_0_x_y_y_yy = buffer_2000_pppd[675];

    auto g_yz_0_0_0_x_y_y_yz = buffer_2000_pppd[676];

    auto g_yz_0_0_0_x_y_y_zz = buffer_2000_pppd[677];

    auto g_yz_0_0_0_x_y_z_xx = buffer_2000_pppd[678];

    auto g_yz_0_0_0_x_y_z_xy = buffer_2000_pppd[679];

    auto g_yz_0_0_0_x_y_z_xz = buffer_2000_pppd[680];

    auto g_yz_0_0_0_x_y_z_yy = buffer_2000_pppd[681];

    auto g_yz_0_0_0_x_y_z_yz = buffer_2000_pppd[682];

    auto g_yz_0_0_0_x_y_z_zz = buffer_2000_pppd[683];

    auto g_yz_0_0_0_x_z_x_xx = buffer_2000_pppd[684];

    auto g_yz_0_0_0_x_z_x_xy = buffer_2000_pppd[685];

    auto g_yz_0_0_0_x_z_x_xz = buffer_2000_pppd[686];

    auto g_yz_0_0_0_x_z_x_yy = buffer_2000_pppd[687];

    auto g_yz_0_0_0_x_z_x_yz = buffer_2000_pppd[688];

    auto g_yz_0_0_0_x_z_x_zz = buffer_2000_pppd[689];

    auto g_yz_0_0_0_x_z_y_xx = buffer_2000_pppd[690];

    auto g_yz_0_0_0_x_z_y_xy = buffer_2000_pppd[691];

    auto g_yz_0_0_0_x_z_y_xz = buffer_2000_pppd[692];

    auto g_yz_0_0_0_x_z_y_yy = buffer_2000_pppd[693];

    auto g_yz_0_0_0_x_z_y_yz = buffer_2000_pppd[694];

    auto g_yz_0_0_0_x_z_y_zz = buffer_2000_pppd[695];

    auto g_yz_0_0_0_x_z_z_xx = buffer_2000_pppd[696];

    auto g_yz_0_0_0_x_z_z_xy = buffer_2000_pppd[697];

    auto g_yz_0_0_0_x_z_z_xz = buffer_2000_pppd[698];

    auto g_yz_0_0_0_x_z_z_yy = buffer_2000_pppd[699];

    auto g_yz_0_0_0_x_z_z_yz = buffer_2000_pppd[700];

    auto g_yz_0_0_0_x_z_z_zz = buffer_2000_pppd[701];

    auto g_yz_0_0_0_y_x_x_xx = buffer_2000_pppd[702];

    auto g_yz_0_0_0_y_x_x_xy = buffer_2000_pppd[703];

    auto g_yz_0_0_0_y_x_x_xz = buffer_2000_pppd[704];

    auto g_yz_0_0_0_y_x_x_yy = buffer_2000_pppd[705];

    auto g_yz_0_0_0_y_x_x_yz = buffer_2000_pppd[706];

    auto g_yz_0_0_0_y_x_x_zz = buffer_2000_pppd[707];

    auto g_yz_0_0_0_y_x_y_xx = buffer_2000_pppd[708];

    auto g_yz_0_0_0_y_x_y_xy = buffer_2000_pppd[709];

    auto g_yz_0_0_0_y_x_y_xz = buffer_2000_pppd[710];

    auto g_yz_0_0_0_y_x_y_yy = buffer_2000_pppd[711];

    auto g_yz_0_0_0_y_x_y_yz = buffer_2000_pppd[712];

    auto g_yz_0_0_0_y_x_y_zz = buffer_2000_pppd[713];

    auto g_yz_0_0_0_y_x_z_xx = buffer_2000_pppd[714];

    auto g_yz_0_0_0_y_x_z_xy = buffer_2000_pppd[715];

    auto g_yz_0_0_0_y_x_z_xz = buffer_2000_pppd[716];

    auto g_yz_0_0_0_y_x_z_yy = buffer_2000_pppd[717];

    auto g_yz_0_0_0_y_x_z_yz = buffer_2000_pppd[718];

    auto g_yz_0_0_0_y_x_z_zz = buffer_2000_pppd[719];

    auto g_yz_0_0_0_y_y_x_xx = buffer_2000_pppd[720];

    auto g_yz_0_0_0_y_y_x_xy = buffer_2000_pppd[721];

    auto g_yz_0_0_0_y_y_x_xz = buffer_2000_pppd[722];

    auto g_yz_0_0_0_y_y_x_yy = buffer_2000_pppd[723];

    auto g_yz_0_0_0_y_y_x_yz = buffer_2000_pppd[724];

    auto g_yz_0_0_0_y_y_x_zz = buffer_2000_pppd[725];

    auto g_yz_0_0_0_y_y_y_xx = buffer_2000_pppd[726];

    auto g_yz_0_0_0_y_y_y_xy = buffer_2000_pppd[727];

    auto g_yz_0_0_0_y_y_y_xz = buffer_2000_pppd[728];

    auto g_yz_0_0_0_y_y_y_yy = buffer_2000_pppd[729];

    auto g_yz_0_0_0_y_y_y_yz = buffer_2000_pppd[730];

    auto g_yz_0_0_0_y_y_y_zz = buffer_2000_pppd[731];

    auto g_yz_0_0_0_y_y_z_xx = buffer_2000_pppd[732];

    auto g_yz_0_0_0_y_y_z_xy = buffer_2000_pppd[733];

    auto g_yz_0_0_0_y_y_z_xz = buffer_2000_pppd[734];

    auto g_yz_0_0_0_y_y_z_yy = buffer_2000_pppd[735];

    auto g_yz_0_0_0_y_y_z_yz = buffer_2000_pppd[736];

    auto g_yz_0_0_0_y_y_z_zz = buffer_2000_pppd[737];

    auto g_yz_0_0_0_y_z_x_xx = buffer_2000_pppd[738];

    auto g_yz_0_0_0_y_z_x_xy = buffer_2000_pppd[739];

    auto g_yz_0_0_0_y_z_x_xz = buffer_2000_pppd[740];

    auto g_yz_0_0_0_y_z_x_yy = buffer_2000_pppd[741];

    auto g_yz_0_0_0_y_z_x_yz = buffer_2000_pppd[742];

    auto g_yz_0_0_0_y_z_x_zz = buffer_2000_pppd[743];

    auto g_yz_0_0_0_y_z_y_xx = buffer_2000_pppd[744];

    auto g_yz_0_0_0_y_z_y_xy = buffer_2000_pppd[745];

    auto g_yz_0_0_0_y_z_y_xz = buffer_2000_pppd[746];

    auto g_yz_0_0_0_y_z_y_yy = buffer_2000_pppd[747];

    auto g_yz_0_0_0_y_z_y_yz = buffer_2000_pppd[748];

    auto g_yz_0_0_0_y_z_y_zz = buffer_2000_pppd[749];

    auto g_yz_0_0_0_y_z_z_xx = buffer_2000_pppd[750];

    auto g_yz_0_0_0_y_z_z_xy = buffer_2000_pppd[751];

    auto g_yz_0_0_0_y_z_z_xz = buffer_2000_pppd[752];

    auto g_yz_0_0_0_y_z_z_yy = buffer_2000_pppd[753];

    auto g_yz_0_0_0_y_z_z_yz = buffer_2000_pppd[754];

    auto g_yz_0_0_0_y_z_z_zz = buffer_2000_pppd[755];

    auto g_yz_0_0_0_z_x_x_xx = buffer_2000_pppd[756];

    auto g_yz_0_0_0_z_x_x_xy = buffer_2000_pppd[757];

    auto g_yz_0_0_0_z_x_x_xz = buffer_2000_pppd[758];

    auto g_yz_0_0_0_z_x_x_yy = buffer_2000_pppd[759];

    auto g_yz_0_0_0_z_x_x_yz = buffer_2000_pppd[760];

    auto g_yz_0_0_0_z_x_x_zz = buffer_2000_pppd[761];

    auto g_yz_0_0_0_z_x_y_xx = buffer_2000_pppd[762];

    auto g_yz_0_0_0_z_x_y_xy = buffer_2000_pppd[763];

    auto g_yz_0_0_0_z_x_y_xz = buffer_2000_pppd[764];

    auto g_yz_0_0_0_z_x_y_yy = buffer_2000_pppd[765];

    auto g_yz_0_0_0_z_x_y_yz = buffer_2000_pppd[766];

    auto g_yz_0_0_0_z_x_y_zz = buffer_2000_pppd[767];

    auto g_yz_0_0_0_z_x_z_xx = buffer_2000_pppd[768];

    auto g_yz_0_0_0_z_x_z_xy = buffer_2000_pppd[769];

    auto g_yz_0_0_0_z_x_z_xz = buffer_2000_pppd[770];

    auto g_yz_0_0_0_z_x_z_yy = buffer_2000_pppd[771];

    auto g_yz_0_0_0_z_x_z_yz = buffer_2000_pppd[772];

    auto g_yz_0_0_0_z_x_z_zz = buffer_2000_pppd[773];

    auto g_yz_0_0_0_z_y_x_xx = buffer_2000_pppd[774];

    auto g_yz_0_0_0_z_y_x_xy = buffer_2000_pppd[775];

    auto g_yz_0_0_0_z_y_x_xz = buffer_2000_pppd[776];

    auto g_yz_0_0_0_z_y_x_yy = buffer_2000_pppd[777];

    auto g_yz_0_0_0_z_y_x_yz = buffer_2000_pppd[778];

    auto g_yz_0_0_0_z_y_x_zz = buffer_2000_pppd[779];

    auto g_yz_0_0_0_z_y_y_xx = buffer_2000_pppd[780];

    auto g_yz_0_0_0_z_y_y_xy = buffer_2000_pppd[781];

    auto g_yz_0_0_0_z_y_y_xz = buffer_2000_pppd[782];

    auto g_yz_0_0_0_z_y_y_yy = buffer_2000_pppd[783];

    auto g_yz_0_0_0_z_y_y_yz = buffer_2000_pppd[784];

    auto g_yz_0_0_0_z_y_y_zz = buffer_2000_pppd[785];

    auto g_yz_0_0_0_z_y_z_xx = buffer_2000_pppd[786];

    auto g_yz_0_0_0_z_y_z_xy = buffer_2000_pppd[787];

    auto g_yz_0_0_0_z_y_z_xz = buffer_2000_pppd[788];

    auto g_yz_0_0_0_z_y_z_yy = buffer_2000_pppd[789];

    auto g_yz_0_0_0_z_y_z_yz = buffer_2000_pppd[790];

    auto g_yz_0_0_0_z_y_z_zz = buffer_2000_pppd[791];

    auto g_yz_0_0_0_z_z_x_xx = buffer_2000_pppd[792];

    auto g_yz_0_0_0_z_z_x_xy = buffer_2000_pppd[793];

    auto g_yz_0_0_0_z_z_x_xz = buffer_2000_pppd[794];

    auto g_yz_0_0_0_z_z_x_yy = buffer_2000_pppd[795];

    auto g_yz_0_0_0_z_z_x_yz = buffer_2000_pppd[796];

    auto g_yz_0_0_0_z_z_x_zz = buffer_2000_pppd[797];

    auto g_yz_0_0_0_z_z_y_xx = buffer_2000_pppd[798];

    auto g_yz_0_0_0_z_z_y_xy = buffer_2000_pppd[799];

    auto g_yz_0_0_0_z_z_y_xz = buffer_2000_pppd[800];

    auto g_yz_0_0_0_z_z_y_yy = buffer_2000_pppd[801];

    auto g_yz_0_0_0_z_z_y_yz = buffer_2000_pppd[802];

    auto g_yz_0_0_0_z_z_y_zz = buffer_2000_pppd[803];

    auto g_yz_0_0_0_z_z_z_xx = buffer_2000_pppd[804];

    auto g_yz_0_0_0_z_z_z_xy = buffer_2000_pppd[805];

    auto g_yz_0_0_0_z_z_z_xz = buffer_2000_pppd[806];

    auto g_yz_0_0_0_z_z_z_yy = buffer_2000_pppd[807];

    auto g_yz_0_0_0_z_z_z_yz = buffer_2000_pppd[808];

    auto g_yz_0_0_0_z_z_z_zz = buffer_2000_pppd[809];

    auto g_zz_0_0_0_x_x_x_xx = buffer_2000_pppd[810];

    auto g_zz_0_0_0_x_x_x_xy = buffer_2000_pppd[811];

    auto g_zz_0_0_0_x_x_x_xz = buffer_2000_pppd[812];

    auto g_zz_0_0_0_x_x_x_yy = buffer_2000_pppd[813];

    auto g_zz_0_0_0_x_x_x_yz = buffer_2000_pppd[814];

    auto g_zz_0_0_0_x_x_x_zz = buffer_2000_pppd[815];

    auto g_zz_0_0_0_x_x_y_xx = buffer_2000_pppd[816];

    auto g_zz_0_0_0_x_x_y_xy = buffer_2000_pppd[817];

    auto g_zz_0_0_0_x_x_y_xz = buffer_2000_pppd[818];

    auto g_zz_0_0_0_x_x_y_yy = buffer_2000_pppd[819];

    auto g_zz_0_0_0_x_x_y_yz = buffer_2000_pppd[820];

    auto g_zz_0_0_0_x_x_y_zz = buffer_2000_pppd[821];

    auto g_zz_0_0_0_x_x_z_xx = buffer_2000_pppd[822];

    auto g_zz_0_0_0_x_x_z_xy = buffer_2000_pppd[823];

    auto g_zz_0_0_0_x_x_z_xz = buffer_2000_pppd[824];

    auto g_zz_0_0_0_x_x_z_yy = buffer_2000_pppd[825];

    auto g_zz_0_0_0_x_x_z_yz = buffer_2000_pppd[826];

    auto g_zz_0_0_0_x_x_z_zz = buffer_2000_pppd[827];

    auto g_zz_0_0_0_x_y_x_xx = buffer_2000_pppd[828];

    auto g_zz_0_0_0_x_y_x_xy = buffer_2000_pppd[829];

    auto g_zz_0_0_0_x_y_x_xz = buffer_2000_pppd[830];

    auto g_zz_0_0_0_x_y_x_yy = buffer_2000_pppd[831];

    auto g_zz_0_0_0_x_y_x_yz = buffer_2000_pppd[832];

    auto g_zz_0_0_0_x_y_x_zz = buffer_2000_pppd[833];

    auto g_zz_0_0_0_x_y_y_xx = buffer_2000_pppd[834];

    auto g_zz_0_0_0_x_y_y_xy = buffer_2000_pppd[835];

    auto g_zz_0_0_0_x_y_y_xz = buffer_2000_pppd[836];

    auto g_zz_0_0_0_x_y_y_yy = buffer_2000_pppd[837];

    auto g_zz_0_0_0_x_y_y_yz = buffer_2000_pppd[838];

    auto g_zz_0_0_0_x_y_y_zz = buffer_2000_pppd[839];

    auto g_zz_0_0_0_x_y_z_xx = buffer_2000_pppd[840];

    auto g_zz_0_0_0_x_y_z_xy = buffer_2000_pppd[841];

    auto g_zz_0_0_0_x_y_z_xz = buffer_2000_pppd[842];

    auto g_zz_0_0_0_x_y_z_yy = buffer_2000_pppd[843];

    auto g_zz_0_0_0_x_y_z_yz = buffer_2000_pppd[844];

    auto g_zz_0_0_0_x_y_z_zz = buffer_2000_pppd[845];

    auto g_zz_0_0_0_x_z_x_xx = buffer_2000_pppd[846];

    auto g_zz_0_0_0_x_z_x_xy = buffer_2000_pppd[847];

    auto g_zz_0_0_0_x_z_x_xz = buffer_2000_pppd[848];

    auto g_zz_0_0_0_x_z_x_yy = buffer_2000_pppd[849];

    auto g_zz_0_0_0_x_z_x_yz = buffer_2000_pppd[850];

    auto g_zz_0_0_0_x_z_x_zz = buffer_2000_pppd[851];

    auto g_zz_0_0_0_x_z_y_xx = buffer_2000_pppd[852];

    auto g_zz_0_0_0_x_z_y_xy = buffer_2000_pppd[853];

    auto g_zz_0_0_0_x_z_y_xz = buffer_2000_pppd[854];

    auto g_zz_0_0_0_x_z_y_yy = buffer_2000_pppd[855];

    auto g_zz_0_0_0_x_z_y_yz = buffer_2000_pppd[856];

    auto g_zz_0_0_0_x_z_y_zz = buffer_2000_pppd[857];

    auto g_zz_0_0_0_x_z_z_xx = buffer_2000_pppd[858];

    auto g_zz_0_0_0_x_z_z_xy = buffer_2000_pppd[859];

    auto g_zz_0_0_0_x_z_z_xz = buffer_2000_pppd[860];

    auto g_zz_0_0_0_x_z_z_yy = buffer_2000_pppd[861];

    auto g_zz_0_0_0_x_z_z_yz = buffer_2000_pppd[862];

    auto g_zz_0_0_0_x_z_z_zz = buffer_2000_pppd[863];

    auto g_zz_0_0_0_y_x_x_xx = buffer_2000_pppd[864];

    auto g_zz_0_0_0_y_x_x_xy = buffer_2000_pppd[865];

    auto g_zz_0_0_0_y_x_x_xz = buffer_2000_pppd[866];

    auto g_zz_0_0_0_y_x_x_yy = buffer_2000_pppd[867];

    auto g_zz_0_0_0_y_x_x_yz = buffer_2000_pppd[868];

    auto g_zz_0_0_0_y_x_x_zz = buffer_2000_pppd[869];

    auto g_zz_0_0_0_y_x_y_xx = buffer_2000_pppd[870];

    auto g_zz_0_0_0_y_x_y_xy = buffer_2000_pppd[871];

    auto g_zz_0_0_0_y_x_y_xz = buffer_2000_pppd[872];

    auto g_zz_0_0_0_y_x_y_yy = buffer_2000_pppd[873];

    auto g_zz_0_0_0_y_x_y_yz = buffer_2000_pppd[874];

    auto g_zz_0_0_0_y_x_y_zz = buffer_2000_pppd[875];

    auto g_zz_0_0_0_y_x_z_xx = buffer_2000_pppd[876];

    auto g_zz_0_0_0_y_x_z_xy = buffer_2000_pppd[877];

    auto g_zz_0_0_0_y_x_z_xz = buffer_2000_pppd[878];

    auto g_zz_0_0_0_y_x_z_yy = buffer_2000_pppd[879];

    auto g_zz_0_0_0_y_x_z_yz = buffer_2000_pppd[880];

    auto g_zz_0_0_0_y_x_z_zz = buffer_2000_pppd[881];

    auto g_zz_0_0_0_y_y_x_xx = buffer_2000_pppd[882];

    auto g_zz_0_0_0_y_y_x_xy = buffer_2000_pppd[883];

    auto g_zz_0_0_0_y_y_x_xz = buffer_2000_pppd[884];

    auto g_zz_0_0_0_y_y_x_yy = buffer_2000_pppd[885];

    auto g_zz_0_0_0_y_y_x_yz = buffer_2000_pppd[886];

    auto g_zz_0_0_0_y_y_x_zz = buffer_2000_pppd[887];

    auto g_zz_0_0_0_y_y_y_xx = buffer_2000_pppd[888];

    auto g_zz_0_0_0_y_y_y_xy = buffer_2000_pppd[889];

    auto g_zz_0_0_0_y_y_y_xz = buffer_2000_pppd[890];

    auto g_zz_0_0_0_y_y_y_yy = buffer_2000_pppd[891];

    auto g_zz_0_0_0_y_y_y_yz = buffer_2000_pppd[892];

    auto g_zz_0_0_0_y_y_y_zz = buffer_2000_pppd[893];

    auto g_zz_0_0_0_y_y_z_xx = buffer_2000_pppd[894];

    auto g_zz_0_0_0_y_y_z_xy = buffer_2000_pppd[895];

    auto g_zz_0_0_0_y_y_z_xz = buffer_2000_pppd[896];

    auto g_zz_0_0_0_y_y_z_yy = buffer_2000_pppd[897];

    auto g_zz_0_0_0_y_y_z_yz = buffer_2000_pppd[898];

    auto g_zz_0_0_0_y_y_z_zz = buffer_2000_pppd[899];

    auto g_zz_0_0_0_y_z_x_xx = buffer_2000_pppd[900];

    auto g_zz_0_0_0_y_z_x_xy = buffer_2000_pppd[901];

    auto g_zz_0_0_0_y_z_x_xz = buffer_2000_pppd[902];

    auto g_zz_0_0_0_y_z_x_yy = buffer_2000_pppd[903];

    auto g_zz_0_0_0_y_z_x_yz = buffer_2000_pppd[904];

    auto g_zz_0_0_0_y_z_x_zz = buffer_2000_pppd[905];

    auto g_zz_0_0_0_y_z_y_xx = buffer_2000_pppd[906];

    auto g_zz_0_0_0_y_z_y_xy = buffer_2000_pppd[907];

    auto g_zz_0_0_0_y_z_y_xz = buffer_2000_pppd[908];

    auto g_zz_0_0_0_y_z_y_yy = buffer_2000_pppd[909];

    auto g_zz_0_0_0_y_z_y_yz = buffer_2000_pppd[910];

    auto g_zz_0_0_0_y_z_y_zz = buffer_2000_pppd[911];

    auto g_zz_0_0_0_y_z_z_xx = buffer_2000_pppd[912];

    auto g_zz_0_0_0_y_z_z_xy = buffer_2000_pppd[913];

    auto g_zz_0_0_0_y_z_z_xz = buffer_2000_pppd[914];

    auto g_zz_0_0_0_y_z_z_yy = buffer_2000_pppd[915];

    auto g_zz_0_0_0_y_z_z_yz = buffer_2000_pppd[916];

    auto g_zz_0_0_0_y_z_z_zz = buffer_2000_pppd[917];

    auto g_zz_0_0_0_z_x_x_xx = buffer_2000_pppd[918];

    auto g_zz_0_0_0_z_x_x_xy = buffer_2000_pppd[919];

    auto g_zz_0_0_0_z_x_x_xz = buffer_2000_pppd[920];

    auto g_zz_0_0_0_z_x_x_yy = buffer_2000_pppd[921];

    auto g_zz_0_0_0_z_x_x_yz = buffer_2000_pppd[922];

    auto g_zz_0_0_0_z_x_x_zz = buffer_2000_pppd[923];

    auto g_zz_0_0_0_z_x_y_xx = buffer_2000_pppd[924];

    auto g_zz_0_0_0_z_x_y_xy = buffer_2000_pppd[925];

    auto g_zz_0_0_0_z_x_y_xz = buffer_2000_pppd[926];

    auto g_zz_0_0_0_z_x_y_yy = buffer_2000_pppd[927];

    auto g_zz_0_0_0_z_x_y_yz = buffer_2000_pppd[928];

    auto g_zz_0_0_0_z_x_y_zz = buffer_2000_pppd[929];

    auto g_zz_0_0_0_z_x_z_xx = buffer_2000_pppd[930];

    auto g_zz_0_0_0_z_x_z_xy = buffer_2000_pppd[931];

    auto g_zz_0_0_0_z_x_z_xz = buffer_2000_pppd[932];

    auto g_zz_0_0_0_z_x_z_yy = buffer_2000_pppd[933];

    auto g_zz_0_0_0_z_x_z_yz = buffer_2000_pppd[934];

    auto g_zz_0_0_0_z_x_z_zz = buffer_2000_pppd[935];

    auto g_zz_0_0_0_z_y_x_xx = buffer_2000_pppd[936];

    auto g_zz_0_0_0_z_y_x_xy = buffer_2000_pppd[937];

    auto g_zz_0_0_0_z_y_x_xz = buffer_2000_pppd[938];

    auto g_zz_0_0_0_z_y_x_yy = buffer_2000_pppd[939];

    auto g_zz_0_0_0_z_y_x_yz = buffer_2000_pppd[940];

    auto g_zz_0_0_0_z_y_x_zz = buffer_2000_pppd[941];

    auto g_zz_0_0_0_z_y_y_xx = buffer_2000_pppd[942];

    auto g_zz_0_0_0_z_y_y_xy = buffer_2000_pppd[943];

    auto g_zz_0_0_0_z_y_y_xz = buffer_2000_pppd[944];

    auto g_zz_0_0_0_z_y_y_yy = buffer_2000_pppd[945];

    auto g_zz_0_0_0_z_y_y_yz = buffer_2000_pppd[946];

    auto g_zz_0_0_0_z_y_y_zz = buffer_2000_pppd[947];

    auto g_zz_0_0_0_z_y_z_xx = buffer_2000_pppd[948];

    auto g_zz_0_0_0_z_y_z_xy = buffer_2000_pppd[949];

    auto g_zz_0_0_0_z_y_z_xz = buffer_2000_pppd[950];

    auto g_zz_0_0_0_z_y_z_yy = buffer_2000_pppd[951];

    auto g_zz_0_0_0_z_y_z_yz = buffer_2000_pppd[952];

    auto g_zz_0_0_0_z_y_z_zz = buffer_2000_pppd[953];

    auto g_zz_0_0_0_z_z_x_xx = buffer_2000_pppd[954];

    auto g_zz_0_0_0_z_z_x_xy = buffer_2000_pppd[955];

    auto g_zz_0_0_0_z_z_x_xz = buffer_2000_pppd[956];

    auto g_zz_0_0_0_z_z_x_yy = buffer_2000_pppd[957];

    auto g_zz_0_0_0_z_z_x_yz = buffer_2000_pppd[958];

    auto g_zz_0_0_0_z_z_x_zz = buffer_2000_pppd[959];

    auto g_zz_0_0_0_z_z_y_xx = buffer_2000_pppd[960];

    auto g_zz_0_0_0_z_z_y_xy = buffer_2000_pppd[961];

    auto g_zz_0_0_0_z_z_y_xz = buffer_2000_pppd[962];

    auto g_zz_0_0_0_z_z_y_yy = buffer_2000_pppd[963];

    auto g_zz_0_0_0_z_z_y_yz = buffer_2000_pppd[964];

    auto g_zz_0_0_0_z_z_y_zz = buffer_2000_pppd[965];

    auto g_zz_0_0_0_z_z_z_xx = buffer_2000_pppd[966];

    auto g_zz_0_0_0_z_z_z_xy = buffer_2000_pppd[967];

    auto g_zz_0_0_0_z_z_z_xz = buffer_2000_pppd[968];

    auto g_zz_0_0_0_z_z_z_yy = buffer_2000_pppd[969];

    auto g_zz_0_0_0_z_z_z_yz = buffer_2000_pppd[970];

    auto g_zz_0_0_0_z_z_z_zz = buffer_2000_pppd[971];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_xx_0_0_0_x_x_x_xx, g_xx_0_0_0_x_x_x_xy, g_xx_0_0_0_x_x_x_xz, g_xx_0_0_0_x_x_x_yy, g_xx_0_0_0_x_x_x_yz, g_xx_0_0_0_x_x_x_zz, g_xxx_x_x_xx, g_xxx_x_x_xy, g_xxx_x_x_xz, g_xxx_x_x_yy, g_xxx_x_x_yz, g_xxx_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_x_xx[i] = -6.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_xxx_x_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_x_xy[i] = -6.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_xxx_x_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_x_xz[i] = -6.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_xxx_x_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_x_yy[i] = -6.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_xxx_x_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_x_yz[i] = -6.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_xxx_x_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_x_zz[i] = -6.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_xxx_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, g_xx_0_0_0_x_x_y_xx, g_xx_0_0_0_x_x_y_xy, g_xx_0_0_0_x_x_y_xz, g_xx_0_0_0_x_x_y_yy, g_xx_0_0_0_x_x_y_yz, g_xx_0_0_0_x_x_y_zz, g_xxx_x_y_xx, g_xxx_x_y_xy, g_xxx_x_y_xz, g_xxx_x_y_yy, g_xxx_x_y_yz, g_xxx_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_y_xx[i] = -6.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_xxx_x_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_y_xy[i] = -6.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_xxx_x_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_y_xz[i] = -6.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_xxx_x_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_y_yy[i] = -6.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_xxx_x_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_y_yz[i] = -6.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_xxx_x_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_y_zz[i] = -6.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_xxx_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, g_xx_0_0_0_x_x_z_xx, g_xx_0_0_0_x_x_z_xy, g_xx_0_0_0_x_x_z_xz, g_xx_0_0_0_x_x_z_yy, g_xx_0_0_0_x_x_z_yz, g_xx_0_0_0_x_x_z_zz, g_xxx_x_z_xx, g_xxx_x_z_xy, g_xxx_x_z_xz, g_xxx_x_z_yy, g_xxx_x_z_yz, g_xxx_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_z_xx[i] = -6.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_xxx_x_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_z_xy[i] = -6.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_xxx_x_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_z_xz[i] = -6.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_xxx_x_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_z_yy[i] = -6.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_xxx_x_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_z_yz[i] = -6.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_xxx_x_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_z_zz[i] = -6.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_xxx_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_xx_0_0_0_x_y_x_xx, g_xx_0_0_0_x_y_x_xy, g_xx_0_0_0_x_y_x_xz, g_xx_0_0_0_x_y_x_yy, g_xx_0_0_0_x_y_x_yz, g_xx_0_0_0_x_y_x_zz, g_xxx_y_x_xx, g_xxx_y_x_xy, g_xxx_y_x_xz, g_xxx_y_x_yy, g_xxx_y_x_yz, g_xxx_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_x_xx[i] = -6.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_xxx_y_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_x_xy[i] = -6.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_xxx_y_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_x_xz[i] = -6.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_xxx_y_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_x_yy[i] = -6.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_xxx_y_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_x_yz[i] = -6.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_xxx_y_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_x_zz[i] = -6.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_xxx_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, g_xx_0_0_0_x_y_y_xx, g_xx_0_0_0_x_y_y_xy, g_xx_0_0_0_x_y_y_xz, g_xx_0_0_0_x_y_y_yy, g_xx_0_0_0_x_y_y_yz, g_xx_0_0_0_x_y_y_zz, g_xxx_y_y_xx, g_xxx_y_y_xy, g_xxx_y_y_xz, g_xxx_y_y_yy, g_xxx_y_y_yz, g_xxx_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_y_xx[i] = -6.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_xxx_y_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_y_xy[i] = -6.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_xxx_y_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_y_xz[i] = -6.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_xxx_y_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_y_yy[i] = -6.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_xxx_y_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_y_yz[i] = -6.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_xxx_y_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_y_zz[i] = -6.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_xxx_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, g_xx_0_0_0_x_y_z_xx, g_xx_0_0_0_x_y_z_xy, g_xx_0_0_0_x_y_z_xz, g_xx_0_0_0_x_y_z_yy, g_xx_0_0_0_x_y_z_yz, g_xx_0_0_0_x_y_z_zz, g_xxx_y_z_xx, g_xxx_y_z_xy, g_xxx_y_z_xz, g_xxx_y_z_yy, g_xxx_y_z_yz, g_xxx_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_z_xx[i] = -6.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_xxx_y_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_z_xy[i] = -6.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_xxx_y_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_z_xz[i] = -6.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_xxx_y_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_z_yy[i] = -6.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_xxx_y_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_z_yz[i] = -6.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_xxx_y_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_z_zz[i] = -6.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_xxx_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, g_xx_0_0_0_x_z_x_xx, g_xx_0_0_0_x_z_x_xy, g_xx_0_0_0_x_z_x_xz, g_xx_0_0_0_x_z_x_yy, g_xx_0_0_0_x_z_x_yz, g_xx_0_0_0_x_z_x_zz, g_xxx_z_x_xx, g_xxx_z_x_xy, g_xxx_z_x_xz, g_xxx_z_x_yy, g_xxx_z_x_yz, g_xxx_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_x_xx[i] = -6.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_xxx_z_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_x_xy[i] = -6.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_xxx_z_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_x_xz[i] = -6.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_xxx_z_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_x_yy[i] = -6.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_xxx_z_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_x_yz[i] = -6.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_xxx_z_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_x_zz[i] = -6.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_xxx_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz, g_xx_0_0_0_x_z_y_xx, g_xx_0_0_0_x_z_y_xy, g_xx_0_0_0_x_z_y_xz, g_xx_0_0_0_x_z_y_yy, g_xx_0_0_0_x_z_y_yz, g_xx_0_0_0_x_z_y_zz, g_xxx_z_y_xx, g_xxx_z_y_xy, g_xxx_z_y_xz, g_xxx_z_y_yy, g_xxx_z_y_yz, g_xxx_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_y_xx[i] = -6.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_xxx_z_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_y_xy[i] = -6.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_xxx_z_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_y_xz[i] = -6.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_xxx_z_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_y_yy[i] = -6.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_xxx_z_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_y_yz[i] = -6.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_xxx_z_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_y_zz[i] = -6.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_xxx_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz, g_xx_0_0_0_x_z_z_xx, g_xx_0_0_0_x_z_z_xy, g_xx_0_0_0_x_z_z_xz, g_xx_0_0_0_x_z_z_yy, g_xx_0_0_0_x_z_z_yz, g_xx_0_0_0_x_z_z_zz, g_xxx_z_z_xx, g_xxx_z_z_xy, g_xxx_z_z_xz, g_xxx_z_z_yy, g_xxx_z_z_yz, g_xxx_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_z_xx[i] = -6.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_xxx_z_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_z_xy[i] = -6.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_xxx_z_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_z_xz[i] = -6.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_xxx_z_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_z_yy[i] = -6.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_xxx_z_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_z_yz[i] = -6.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_xxx_z_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_z_zz[i] = -6.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_xxx_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_x_xx, g_xx_0_0_0_y_x_x_xy, g_xx_0_0_0_y_x_x_xz, g_xx_0_0_0_y_x_x_yy, g_xx_0_0_0_y_x_x_yz, g_xx_0_0_0_y_x_x_zz, g_xxy_x_x_xx, g_xxy_x_x_xy, g_xxy_x_x_xz, g_xxy_x_x_yy, g_xxy_x_x_yz, g_xxy_x_x_zz, g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_x_xx[i] = -2.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_xxy_x_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_x_xy[i] = -2.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_xxy_x_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_x_xz[i] = -2.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_xxy_x_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_x_yy[i] = -2.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_xxy_x_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_x_yz[i] = -2.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_xxy_x_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_x_zz[i] = -2.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_xxy_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_y_xx, g_xx_0_0_0_y_x_y_xy, g_xx_0_0_0_y_x_y_xz, g_xx_0_0_0_y_x_y_yy, g_xx_0_0_0_y_x_y_yz, g_xx_0_0_0_y_x_y_zz, g_xxy_x_y_xx, g_xxy_x_y_xy, g_xxy_x_y_xz, g_xxy_x_y_yy, g_xxy_x_y_yz, g_xxy_x_y_zz, g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_y_xx[i] = -2.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_xxy_x_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_y_xy[i] = -2.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_xxy_x_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_y_xz[i] = -2.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_xxy_x_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_y_yy[i] = -2.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_xxy_x_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_y_yz[i] = -2.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_xxy_x_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_y_zz[i] = -2.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_xxy_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_z_xx, g_xx_0_0_0_y_x_z_xy, g_xx_0_0_0_y_x_z_xz, g_xx_0_0_0_y_x_z_yy, g_xx_0_0_0_y_x_z_yz, g_xx_0_0_0_y_x_z_zz, g_xxy_x_z_xx, g_xxy_x_z_xy, g_xxy_x_z_xz, g_xxy_x_z_yy, g_xxy_x_z_yz, g_xxy_x_z_zz, g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_z_xx[i] = -2.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_xxy_x_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_z_xy[i] = -2.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_xxy_x_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_z_xz[i] = -2.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_xxy_x_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_z_yy[i] = -2.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_xxy_x_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_z_yz[i] = -2.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_xxy_x_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_z_zz[i] = -2.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_xxy_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_x_xx, g_xx_0_0_0_y_y_x_xy, g_xx_0_0_0_y_y_x_xz, g_xx_0_0_0_y_y_x_yy, g_xx_0_0_0_y_y_x_yz, g_xx_0_0_0_y_y_x_zz, g_xxy_y_x_xx, g_xxy_y_x_xy, g_xxy_y_x_xz, g_xxy_y_x_yy, g_xxy_y_x_yz, g_xxy_y_x_zz, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_x_xx[i] = -2.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_xxy_y_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_x_xy[i] = -2.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_xxy_y_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_x_xz[i] = -2.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_xxy_y_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_x_yy[i] = -2.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_xxy_y_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_x_yz[i] = -2.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_xxy_y_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_x_zz[i] = -2.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_xxy_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_y_xx, g_xx_0_0_0_y_y_y_xy, g_xx_0_0_0_y_y_y_xz, g_xx_0_0_0_y_y_y_yy, g_xx_0_0_0_y_y_y_yz, g_xx_0_0_0_y_y_y_zz, g_xxy_y_y_xx, g_xxy_y_y_xy, g_xxy_y_y_xz, g_xxy_y_y_yy, g_xxy_y_y_yz, g_xxy_y_y_zz, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_y_xx[i] = -2.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_xxy_y_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_y_xy[i] = -2.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_xxy_y_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_y_xz[i] = -2.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_xxy_y_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_y_yy[i] = -2.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_xxy_y_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_y_yz[i] = -2.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_xxy_y_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_y_zz[i] = -2.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_xxy_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_z_xx, g_xx_0_0_0_y_y_z_xy, g_xx_0_0_0_y_y_z_xz, g_xx_0_0_0_y_y_z_yy, g_xx_0_0_0_y_y_z_yz, g_xx_0_0_0_y_y_z_zz, g_xxy_y_z_xx, g_xxy_y_z_xy, g_xxy_y_z_xz, g_xxy_y_z_yy, g_xxy_y_z_yz, g_xxy_y_z_zz, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_z_xx[i] = -2.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_xxy_y_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_z_xy[i] = -2.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_xxy_y_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_z_xz[i] = -2.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_xxy_y_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_z_yy[i] = -2.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_xxy_y_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_z_yz[i] = -2.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_xxy_y_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_z_zz[i] = -2.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_xxy_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_x_xx, g_xx_0_0_0_y_z_x_xy, g_xx_0_0_0_y_z_x_xz, g_xx_0_0_0_y_z_x_yy, g_xx_0_0_0_y_z_x_yz, g_xx_0_0_0_y_z_x_zz, g_xxy_z_x_xx, g_xxy_z_x_xy, g_xxy_z_x_xz, g_xxy_z_x_yy, g_xxy_z_x_yz, g_xxy_z_x_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_x_xx[i] = -2.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_xxy_z_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_x_xy[i] = -2.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_xxy_z_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_x_xz[i] = -2.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_xxy_z_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_x_yy[i] = -2.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_xxy_z_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_x_yz[i] = -2.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_xxy_z_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_x_zz[i] = -2.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_xxy_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_y_xx, g_xx_0_0_0_y_z_y_xy, g_xx_0_0_0_y_z_y_xz, g_xx_0_0_0_y_z_y_yy, g_xx_0_0_0_y_z_y_yz, g_xx_0_0_0_y_z_y_zz, g_xxy_z_y_xx, g_xxy_z_y_xy, g_xxy_z_y_xz, g_xxy_z_y_yy, g_xxy_z_y_yz, g_xxy_z_y_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_y_xx[i] = -2.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_xxy_z_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_y_xy[i] = -2.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_xxy_z_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_y_xz[i] = -2.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_xxy_z_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_y_yy[i] = -2.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_xxy_z_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_y_yz[i] = -2.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_xxy_z_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_y_zz[i] = -2.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_xxy_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_z_xx, g_xx_0_0_0_y_z_z_xy, g_xx_0_0_0_y_z_z_xz, g_xx_0_0_0_y_z_z_yy, g_xx_0_0_0_y_z_z_yz, g_xx_0_0_0_y_z_z_zz, g_xxy_z_z_xx, g_xxy_z_z_xy, g_xxy_z_z_xz, g_xxy_z_z_yy, g_xxy_z_z_yz, g_xxy_z_z_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_z_xx[i] = -2.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_xxy_z_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_z_xy[i] = -2.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_xxy_z_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_z_xz[i] = -2.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_xxy_z_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_z_yy[i] = -2.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_xxy_z_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_z_yz[i] = -2.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_xxy_z_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_z_zz[i] = -2.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_xxy_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_x_xx, g_xx_0_0_0_z_x_x_xy, g_xx_0_0_0_z_x_x_xz, g_xx_0_0_0_z_x_x_yy, g_xx_0_0_0_z_x_x_yz, g_xx_0_0_0_z_x_x_zz, g_xxz_x_x_xx, g_xxz_x_x_xy, g_xxz_x_x_xz, g_xxz_x_x_yy, g_xxz_x_x_yz, g_xxz_x_x_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_x_xx[i] = -2.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_xxz_x_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_x_xy[i] = -2.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_xxz_x_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_x_xz[i] = -2.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_xxz_x_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_x_yy[i] = -2.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_xxz_x_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_x_yz[i] = -2.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_xxz_x_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_x_zz[i] = -2.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_xxz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_y_xx, g_xx_0_0_0_z_x_y_xy, g_xx_0_0_0_z_x_y_xz, g_xx_0_0_0_z_x_y_yy, g_xx_0_0_0_z_x_y_yz, g_xx_0_0_0_z_x_y_zz, g_xxz_x_y_xx, g_xxz_x_y_xy, g_xxz_x_y_xz, g_xxz_x_y_yy, g_xxz_x_y_yz, g_xxz_x_y_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_y_xx[i] = -2.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_xxz_x_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_y_xy[i] = -2.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_xxz_x_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_y_xz[i] = -2.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_xxz_x_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_y_yy[i] = -2.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_xxz_x_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_y_yz[i] = -2.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_xxz_x_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_y_zz[i] = -2.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_xxz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_z_xx, g_xx_0_0_0_z_x_z_xy, g_xx_0_0_0_z_x_z_xz, g_xx_0_0_0_z_x_z_yy, g_xx_0_0_0_z_x_z_yz, g_xx_0_0_0_z_x_z_zz, g_xxz_x_z_xx, g_xxz_x_z_xy, g_xxz_x_z_xz, g_xxz_x_z_yy, g_xxz_x_z_yz, g_xxz_x_z_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_z_xx[i] = -2.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_xxz_x_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_z_xy[i] = -2.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_xxz_x_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_z_xz[i] = -2.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_xxz_x_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_z_yy[i] = -2.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_xxz_x_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_z_yz[i] = -2.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_xxz_x_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_z_zz[i] = -2.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_xxz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_x_xx, g_xx_0_0_0_z_y_x_xy, g_xx_0_0_0_z_y_x_xz, g_xx_0_0_0_z_y_x_yy, g_xx_0_0_0_z_y_x_yz, g_xx_0_0_0_z_y_x_zz, g_xxz_y_x_xx, g_xxz_y_x_xy, g_xxz_y_x_xz, g_xxz_y_x_yy, g_xxz_y_x_yz, g_xxz_y_x_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_x_xx[i] = -2.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_xxz_y_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_x_xy[i] = -2.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_xxz_y_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_x_xz[i] = -2.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_xxz_y_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_x_yy[i] = -2.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_xxz_y_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_x_yz[i] = -2.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_xxz_y_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_x_zz[i] = -2.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_xxz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_y_xx, g_xx_0_0_0_z_y_y_xy, g_xx_0_0_0_z_y_y_xz, g_xx_0_0_0_z_y_y_yy, g_xx_0_0_0_z_y_y_yz, g_xx_0_0_0_z_y_y_zz, g_xxz_y_y_xx, g_xxz_y_y_xy, g_xxz_y_y_xz, g_xxz_y_y_yy, g_xxz_y_y_yz, g_xxz_y_y_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_y_xx[i] = -2.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_xxz_y_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_y_xy[i] = -2.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_xxz_y_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_y_xz[i] = -2.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_xxz_y_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_y_yy[i] = -2.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_xxz_y_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_y_yz[i] = -2.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_xxz_y_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_y_zz[i] = -2.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_xxz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_z_xx, g_xx_0_0_0_z_y_z_xy, g_xx_0_0_0_z_y_z_xz, g_xx_0_0_0_z_y_z_yy, g_xx_0_0_0_z_y_z_yz, g_xx_0_0_0_z_y_z_zz, g_xxz_y_z_xx, g_xxz_y_z_xy, g_xxz_y_z_xz, g_xxz_y_z_yy, g_xxz_y_z_yz, g_xxz_y_z_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_z_xx[i] = -2.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_xxz_y_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_z_xy[i] = -2.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_xxz_y_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_z_xz[i] = -2.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_xxz_y_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_z_yy[i] = -2.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_xxz_y_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_z_yz[i] = -2.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_xxz_y_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_z_zz[i] = -2.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_xxz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_x_xx, g_xx_0_0_0_z_z_x_xy, g_xx_0_0_0_z_z_x_xz, g_xx_0_0_0_z_z_x_yy, g_xx_0_0_0_z_z_x_yz, g_xx_0_0_0_z_z_x_zz, g_xxz_z_x_xx, g_xxz_z_x_xy, g_xxz_z_x_xz, g_xxz_z_x_yy, g_xxz_z_x_yz, g_xxz_z_x_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_x_xx[i] = -2.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_xxz_z_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_x_xy[i] = -2.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_xxz_z_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_x_xz[i] = -2.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_xxz_z_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_x_yy[i] = -2.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_xxz_z_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_x_yz[i] = -2.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_xxz_z_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_x_zz[i] = -2.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_xxz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_y_xx, g_xx_0_0_0_z_z_y_xy, g_xx_0_0_0_z_z_y_xz, g_xx_0_0_0_z_z_y_yy, g_xx_0_0_0_z_z_y_yz, g_xx_0_0_0_z_z_y_zz, g_xxz_z_y_xx, g_xxz_z_y_xy, g_xxz_z_y_xz, g_xxz_z_y_yy, g_xxz_z_y_yz, g_xxz_z_y_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_y_xx[i] = -2.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_xxz_z_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_y_xy[i] = -2.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_xxz_z_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_y_xz[i] = -2.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_xxz_z_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_y_yy[i] = -2.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_xxz_z_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_y_yz[i] = -2.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_xxz_z_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_y_zz[i] = -2.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_xxz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_z_xx, g_xx_0_0_0_z_z_z_xy, g_xx_0_0_0_z_z_z_xz, g_xx_0_0_0_z_z_z_yy, g_xx_0_0_0_z_z_z_yz, g_xx_0_0_0_z_z_z_zz, g_xxz_z_z_xx, g_xxz_z_z_xy, g_xxz_z_z_xz, g_xxz_z_z_yy, g_xxz_z_z_yz, g_xxz_z_z_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_z_xx[i] = -2.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_xxz_z_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_z_xy[i] = -2.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_xxz_z_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_z_xz[i] = -2.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_xxz_z_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_z_yy[i] = -2.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_xxz_z_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_z_yz[i] = -2.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_xxz_z_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_z_zz[i] = -2.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_xxz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_xxy_x_x_xx, g_xxy_x_x_xy, g_xxy_x_x_xz, g_xxy_x_x_yy, g_xxy_x_x_yz, g_xxy_x_x_zz, g_xy_0_0_0_x_x_x_xx, g_xy_0_0_0_x_x_x_xy, g_xy_0_0_0_x_x_x_xz, g_xy_0_0_0_x_x_x_yy, g_xy_0_0_0_x_x_x_yz, g_xy_0_0_0_x_x_x_zz, g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_x_xx[i] = -2.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_xxy_x_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_x_xy[i] = -2.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_xxy_x_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_x_xz[i] = -2.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_xxy_x_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_x_yy[i] = -2.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_xxy_x_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_x_yz[i] = -2.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_xxy_x_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_x_zz[i] = -2.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_xxy_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_xxy_x_y_xx, g_xxy_x_y_xy, g_xxy_x_y_xz, g_xxy_x_y_yy, g_xxy_x_y_yz, g_xxy_x_y_zz, g_xy_0_0_0_x_x_y_xx, g_xy_0_0_0_x_x_y_xy, g_xy_0_0_0_x_x_y_xz, g_xy_0_0_0_x_x_y_yy, g_xy_0_0_0_x_x_y_yz, g_xy_0_0_0_x_x_y_zz, g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_y_xx[i] = -2.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_xxy_x_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_y_xy[i] = -2.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_xxy_x_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_y_xz[i] = -2.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_xxy_x_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_y_yy[i] = -2.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_xxy_x_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_y_yz[i] = -2.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_xxy_x_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_y_zz[i] = -2.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_xxy_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_xxy_x_z_xx, g_xxy_x_z_xy, g_xxy_x_z_xz, g_xxy_x_z_yy, g_xxy_x_z_yz, g_xxy_x_z_zz, g_xy_0_0_0_x_x_z_xx, g_xy_0_0_0_x_x_z_xy, g_xy_0_0_0_x_x_z_xz, g_xy_0_0_0_x_x_z_yy, g_xy_0_0_0_x_x_z_yz, g_xy_0_0_0_x_x_z_zz, g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_z_xx[i] = -2.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_xxy_x_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_z_xy[i] = -2.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_xxy_x_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_z_xz[i] = -2.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_xxy_x_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_z_yy[i] = -2.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_xxy_x_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_z_yz[i] = -2.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_xxy_x_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_z_zz[i] = -2.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_xxy_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xxy_y_x_xx, g_xxy_y_x_xy, g_xxy_y_x_xz, g_xxy_y_x_yy, g_xxy_y_x_yz, g_xxy_y_x_zz, g_xy_0_0_0_x_y_x_xx, g_xy_0_0_0_x_y_x_xy, g_xy_0_0_0_x_y_x_xz, g_xy_0_0_0_x_y_x_yy, g_xy_0_0_0_x_y_x_yz, g_xy_0_0_0_x_y_x_zz, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_x_xx[i] = -2.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_xxy_y_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_x_xy[i] = -2.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_xxy_y_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_x_xz[i] = -2.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_xxy_y_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_x_yy[i] = -2.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_xxy_y_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_x_yz[i] = -2.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_xxy_y_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_x_zz[i] = -2.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_xxy_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_xxy_y_y_xx, g_xxy_y_y_xy, g_xxy_y_y_xz, g_xxy_y_y_yy, g_xxy_y_y_yz, g_xxy_y_y_zz, g_xy_0_0_0_x_y_y_xx, g_xy_0_0_0_x_y_y_xy, g_xy_0_0_0_x_y_y_xz, g_xy_0_0_0_x_y_y_yy, g_xy_0_0_0_x_y_y_yz, g_xy_0_0_0_x_y_y_zz, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_y_xx[i] = -2.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_xxy_y_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_y_xy[i] = -2.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_xxy_y_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_y_xz[i] = -2.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_xxy_y_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_y_yy[i] = -2.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_xxy_y_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_y_yz[i] = -2.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_xxy_y_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_y_zz[i] = -2.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_xxy_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xxy_y_z_xx, g_xxy_y_z_xy, g_xxy_y_z_xz, g_xxy_y_z_yy, g_xxy_y_z_yz, g_xxy_y_z_zz, g_xy_0_0_0_x_y_z_xx, g_xy_0_0_0_x_y_z_xy, g_xy_0_0_0_x_y_z_xz, g_xy_0_0_0_x_y_z_yy, g_xy_0_0_0_x_y_z_yz, g_xy_0_0_0_x_y_z_zz, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_z_xx[i] = -2.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_xxy_y_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_z_xy[i] = -2.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_xxy_y_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_z_xz[i] = -2.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_xxy_y_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_z_yy[i] = -2.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_xxy_y_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_z_yz[i] = -2.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_xxy_y_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_z_zz[i] = -2.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_xxy_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_xxy_z_x_xx, g_xxy_z_x_xy, g_xxy_z_x_xz, g_xxy_z_x_yy, g_xxy_z_x_yz, g_xxy_z_x_zz, g_xy_0_0_0_x_z_x_xx, g_xy_0_0_0_x_z_x_xy, g_xy_0_0_0_x_z_x_xz, g_xy_0_0_0_x_z_x_yy, g_xy_0_0_0_x_z_x_yz, g_xy_0_0_0_x_z_x_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_x_xx[i] = -2.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_xxy_z_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_x_xy[i] = -2.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_xxy_z_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_x_xz[i] = -2.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_xxy_z_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_x_yy[i] = -2.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_xxy_z_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_x_yz[i] = -2.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_xxy_z_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_x_zz[i] = -2.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_xxy_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_xxy_z_y_xx, g_xxy_z_y_xy, g_xxy_z_y_xz, g_xxy_z_y_yy, g_xxy_z_y_yz, g_xxy_z_y_zz, g_xy_0_0_0_x_z_y_xx, g_xy_0_0_0_x_z_y_xy, g_xy_0_0_0_x_z_y_xz, g_xy_0_0_0_x_z_y_yy, g_xy_0_0_0_x_z_y_yz, g_xy_0_0_0_x_z_y_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_y_xx[i] = -2.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_xxy_z_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_y_xy[i] = -2.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_xxy_z_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_y_xz[i] = -2.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_xxy_z_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_y_yy[i] = -2.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_xxy_z_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_y_yz[i] = -2.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_xxy_z_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_y_zz[i] = -2.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_xxy_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_xxy_z_z_xx, g_xxy_z_z_xy, g_xxy_z_z_xz, g_xxy_z_z_yy, g_xxy_z_z_yz, g_xxy_z_z_zz, g_xy_0_0_0_x_z_z_xx, g_xy_0_0_0_x_z_z_xy, g_xy_0_0_0_x_z_z_xz, g_xy_0_0_0_x_z_z_yy, g_xy_0_0_0_x_z_z_yz, g_xy_0_0_0_x_z_z_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_z_xx[i] = -2.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_xxy_z_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_z_xy[i] = -2.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_xxy_z_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_z_xz[i] = -2.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_xxy_z_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_z_yy[i] = -2.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_xxy_z_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_z_yz[i] = -2.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_xxy_z_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_z_zz[i] = -2.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_xxy_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_xy_0_0_0_y_x_x_xx, g_xy_0_0_0_y_x_x_xy, g_xy_0_0_0_y_x_x_xz, g_xy_0_0_0_y_x_x_yy, g_xy_0_0_0_y_x_x_yz, g_xy_0_0_0_y_x_x_zz, g_xyy_x_x_xx, g_xyy_x_x_xy, g_xyy_x_x_xz, g_xyy_x_x_yy, g_xyy_x_x_yz, g_xyy_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_x_xx[i] = -2.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_xyy_x_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_x_xy[i] = -2.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_xyy_x_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_x_xz[i] = -2.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_xyy_x_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_x_yy[i] = -2.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_xyy_x_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_x_yz[i] = -2.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_xyy_x_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_x_zz[i] = -2.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_xyy_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, g_xy_0_0_0_y_x_y_xx, g_xy_0_0_0_y_x_y_xy, g_xy_0_0_0_y_x_y_xz, g_xy_0_0_0_y_x_y_yy, g_xy_0_0_0_y_x_y_yz, g_xy_0_0_0_y_x_y_zz, g_xyy_x_y_xx, g_xyy_x_y_xy, g_xyy_x_y_xz, g_xyy_x_y_yy, g_xyy_x_y_yz, g_xyy_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_y_xx[i] = -2.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_xyy_x_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_y_xy[i] = -2.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_xyy_x_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_y_xz[i] = -2.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_xyy_x_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_y_yy[i] = -2.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_xyy_x_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_y_yz[i] = -2.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_xyy_x_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_y_zz[i] = -2.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_xyy_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, g_xy_0_0_0_y_x_z_xx, g_xy_0_0_0_y_x_z_xy, g_xy_0_0_0_y_x_z_xz, g_xy_0_0_0_y_x_z_yy, g_xy_0_0_0_y_x_z_yz, g_xy_0_0_0_y_x_z_zz, g_xyy_x_z_xx, g_xyy_x_z_xy, g_xyy_x_z_xz, g_xyy_x_z_yy, g_xyy_x_z_yz, g_xyy_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_z_xx[i] = -2.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_xyy_x_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_z_xy[i] = -2.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_xyy_x_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_z_xz[i] = -2.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_xyy_x_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_z_yy[i] = -2.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_xyy_x_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_z_yz[i] = -2.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_xyy_x_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_z_zz[i] = -2.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_xyy_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_xy_0_0_0_y_y_x_xx, g_xy_0_0_0_y_y_x_xy, g_xy_0_0_0_y_y_x_xz, g_xy_0_0_0_y_y_x_yy, g_xy_0_0_0_y_y_x_yz, g_xy_0_0_0_y_y_x_zz, g_xyy_y_x_xx, g_xyy_y_x_xy, g_xyy_y_x_xz, g_xyy_y_x_yy, g_xyy_y_x_yz, g_xyy_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_x_xx[i] = -2.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_xyy_y_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_x_xy[i] = -2.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_xyy_y_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_x_xz[i] = -2.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_xyy_y_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_x_yy[i] = -2.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_xyy_y_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_x_yz[i] = -2.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_xyy_y_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_x_zz[i] = -2.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_xyy_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, g_xy_0_0_0_y_y_y_xx, g_xy_0_0_0_y_y_y_xy, g_xy_0_0_0_y_y_y_xz, g_xy_0_0_0_y_y_y_yy, g_xy_0_0_0_y_y_y_yz, g_xy_0_0_0_y_y_y_zz, g_xyy_y_y_xx, g_xyy_y_y_xy, g_xyy_y_y_xz, g_xyy_y_y_yy, g_xyy_y_y_yz, g_xyy_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_y_xx[i] = -2.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_xyy_y_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_y_xy[i] = -2.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_xyy_y_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_y_xz[i] = -2.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_xyy_y_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_y_yy[i] = -2.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_xyy_y_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_y_yz[i] = -2.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_xyy_y_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_y_zz[i] = -2.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_xyy_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, g_xy_0_0_0_y_y_z_xx, g_xy_0_0_0_y_y_z_xy, g_xy_0_0_0_y_y_z_xz, g_xy_0_0_0_y_y_z_yy, g_xy_0_0_0_y_y_z_yz, g_xy_0_0_0_y_y_z_zz, g_xyy_y_z_xx, g_xyy_y_z_xy, g_xyy_y_z_xz, g_xyy_y_z_yy, g_xyy_y_z_yz, g_xyy_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_z_xx[i] = -2.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_xyy_y_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_z_xy[i] = -2.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_xyy_y_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_z_xz[i] = -2.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_xyy_y_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_z_yy[i] = -2.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_xyy_y_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_z_yz[i] = -2.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_xyy_y_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_z_zz[i] = -2.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_xyy_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, g_xy_0_0_0_y_z_x_xx, g_xy_0_0_0_y_z_x_xy, g_xy_0_0_0_y_z_x_xz, g_xy_0_0_0_y_z_x_yy, g_xy_0_0_0_y_z_x_yz, g_xy_0_0_0_y_z_x_zz, g_xyy_z_x_xx, g_xyy_z_x_xy, g_xyy_z_x_xz, g_xyy_z_x_yy, g_xyy_z_x_yz, g_xyy_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_x_xx[i] = -2.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_xyy_z_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_x_xy[i] = -2.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_xyy_z_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_x_xz[i] = -2.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_xyy_z_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_x_yy[i] = -2.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_xyy_z_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_x_yz[i] = -2.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_xyy_z_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_x_zz[i] = -2.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_xyy_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz, g_xy_0_0_0_y_z_y_xx, g_xy_0_0_0_y_z_y_xy, g_xy_0_0_0_y_z_y_xz, g_xy_0_0_0_y_z_y_yy, g_xy_0_0_0_y_z_y_yz, g_xy_0_0_0_y_z_y_zz, g_xyy_z_y_xx, g_xyy_z_y_xy, g_xyy_z_y_xz, g_xyy_z_y_yy, g_xyy_z_y_yz, g_xyy_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_y_xx[i] = -2.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_xyy_z_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_y_xy[i] = -2.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_xyy_z_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_y_xz[i] = -2.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_xyy_z_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_y_yy[i] = -2.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_xyy_z_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_y_yz[i] = -2.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_xyy_z_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_y_zz[i] = -2.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_xyy_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz, g_xy_0_0_0_y_z_z_xx, g_xy_0_0_0_y_z_z_xy, g_xy_0_0_0_y_z_z_xz, g_xy_0_0_0_y_z_z_yy, g_xy_0_0_0_y_z_z_yz, g_xy_0_0_0_y_z_z_zz, g_xyy_z_z_xx, g_xyy_z_z_xy, g_xyy_z_z_xz, g_xyy_z_z_yy, g_xyy_z_z_yz, g_xyy_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_z_xx[i] = -2.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_xyy_z_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_z_xy[i] = -2.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_xyy_z_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_z_xz[i] = -2.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_xyy_z_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_z_yy[i] = -2.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_xyy_z_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_z_yz[i] = -2.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_xyy_z_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_z_zz[i] = -2.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_xyy_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_x_xx, g_xy_0_0_0_z_x_x_xy, g_xy_0_0_0_z_x_x_xz, g_xy_0_0_0_z_x_x_yy, g_xy_0_0_0_z_x_x_yz, g_xy_0_0_0_z_x_x_zz, g_xyz_x_x_xx, g_xyz_x_x_xy, g_xyz_x_x_xz, g_xyz_x_x_yy, g_xyz_x_x_yz, g_xyz_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_x_xx[i] = 4.0 * g_xyz_x_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_x_xy[i] = 4.0 * g_xyz_x_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_x_xz[i] = 4.0 * g_xyz_x_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_x_yy[i] = 4.0 * g_xyz_x_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_x_yz[i] = 4.0 * g_xyz_x_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_x_zz[i] = 4.0 * g_xyz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_y_xx, g_xy_0_0_0_z_x_y_xy, g_xy_0_0_0_z_x_y_xz, g_xy_0_0_0_z_x_y_yy, g_xy_0_0_0_z_x_y_yz, g_xy_0_0_0_z_x_y_zz, g_xyz_x_y_xx, g_xyz_x_y_xy, g_xyz_x_y_xz, g_xyz_x_y_yy, g_xyz_x_y_yz, g_xyz_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_y_xx[i] = 4.0 * g_xyz_x_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_y_xy[i] = 4.0 * g_xyz_x_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_y_xz[i] = 4.0 * g_xyz_x_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_y_yy[i] = 4.0 * g_xyz_x_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_y_yz[i] = 4.0 * g_xyz_x_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_y_zz[i] = 4.0 * g_xyz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_z_xx, g_xy_0_0_0_z_x_z_xy, g_xy_0_0_0_z_x_z_xz, g_xy_0_0_0_z_x_z_yy, g_xy_0_0_0_z_x_z_yz, g_xy_0_0_0_z_x_z_zz, g_xyz_x_z_xx, g_xyz_x_z_xy, g_xyz_x_z_xz, g_xyz_x_z_yy, g_xyz_x_z_yz, g_xyz_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_z_xx[i] = 4.0 * g_xyz_x_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_z_xy[i] = 4.0 * g_xyz_x_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_z_xz[i] = 4.0 * g_xyz_x_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_z_yy[i] = 4.0 * g_xyz_x_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_z_yz[i] = 4.0 * g_xyz_x_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_z_zz[i] = 4.0 * g_xyz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_x_xx, g_xy_0_0_0_z_y_x_xy, g_xy_0_0_0_z_y_x_xz, g_xy_0_0_0_z_y_x_yy, g_xy_0_0_0_z_y_x_yz, g_xy_0_0_0_z_y_x_zz, g_xyz_y_x_xx, g_xyz_y_x_xy, g_xyz_y_x_xz, g_xyz_y_x_yy, g_xyz_y_x_yz, g_xyz_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_x_xx[i] = 4.0 * g_xyz_y_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_x_xy[i] = 4.0 * g_xyz_y_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_x_xz[i] = 4.0 * g_xyz_y_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_x_yy[i] = 4.0 * g_xyz_y_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_x_yz[i] = 4.0 * g_xyz_y_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_x_zz[i] = 4.0 * g_xyz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_y_xx, g_xy_0_0_0_z_y_y_xy, g_xy_0_0_0_z_y_y_xz, g_xy_0_0_0_z_y_y_yy, g_xy_0_0_0_z_y_y_yz, g_xy_0_0_0_z_y_y_zz, g_xyz_y_y_xx, g_xyz_y_y_xy, g_xyz_y_y_xz, g_xyz_y_y_yy, g_xyz_y_y_yz, g_xyz_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_y_xx[i] = 4.0 * g_xyz_y_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_y_xy[i] = 4.0 * g_xyz_y_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_y_xz[i] = 4.0 * g_xyz_y_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_y_yy[i] = 4.0 * g_xyz_y_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_y_yz[i] = 4.0 * g_xyz_y_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_y_zz[i] = 4.0 * g_xyz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_z_xx, g_xy_0_0_0_z_y_z_xy, g_xy_0_0_0_z_y_z_xz, g_xy_0_0_0_z_y_z_yy, g_xy_0_0_0_z_y_z_yz, g_xy_0_0_0_z_y_z_zz, g_xyz_y_z_xx, g_xyz_y_z_xy, g_xyz_y_z_xz, g_xyz_y_z_yy, g_xyz_y_z_yz, g_xyz_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_z_xx[i] = 4.0 * g_xyz_y_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_z_xy[i] = 4.0 * g_xyz_y_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_z_xz[i] = 4.0 * g_xyz_y_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_z_yy[i] = 4.0 * g_xyz_y_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_z_yz[i] = 4.0 * g_xyz_y_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_z_zz[i] = 4.0 * g_xyz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_x_xx, g_xy_0_0_0_z_z_x_xy, g_xy_0_0_0_z_z_x_xz, g_xy_0_0_0_z_z_x_yy, g_xy_0_0_0_z_z_x_yz, g_xy_0_0_0_z_z_x_zz, g_xyz_z_x_xx, g_xyz_z_x_xy, g_xyz_z_x_xz, g_xyz_z_x_yy, g_xyz_z_x_yz, g_xyz_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_x_xx[i] = 4.0 * g_xyz_z_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_x_xy[i] = 4.0 * g_xyz_z_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_x_xz[i] = 4.0 * g_xyz_z_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_x_yy[i] = 4.0 * g_xyz_z_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_x_yz[i] = 4.0 * g_xyz_z_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_x_zz[i] = 4.0 * g_xyz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_y_xx, g_xy_0_0_0_z_z_y_xy, g_xy_0_0_0_z_z_y_xz, g_xy_0_0_0_z_z_y_yy, g_xy_0_0_0_z_z_y_yz, g_xy_0_0_0_z_z_y_zz, g_xyz_z_y_xx, g_xyz_z_y_xy, g_xyz_z_y_xz, g_xyz_z_y_yy, g_xyz_z_y_yz, g_xyz_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_y_xx[i] = 4.0 * g_xyz_z_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_y_xy[i] = 4.0 * g_xyz_z_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_y_xz[i] = 4.0 * g_xyz_z_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_y_yy[i] = 4.0 * g_xyz_z_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_y_yz[i] = 4.0 * g_xyz_z_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_y_zz[i] = 4.0 * g_xyz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_z_xx, g_xy_0_0_0_z_z_z_xy, g_xy_0_0_0_z_z_z_xz, g_xy_0_0_0_z_z_z_yy, g_xy_0_0_0_z_z_z_yz, g_xy_0_0_0_z_z_z_zz, g_xyz_z_z_xx, g_xyz_z_z_xy, g_xyz_z_z_xz, g_xyz_z_z_yy, g_xyz_z_z_yz, g_xyz_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_z_xx[i] = 4.0 * g_xyz_z_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_z_xy[i] = 4.0 * g_xyz_z_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_z_xz[i] = 4.0 * g_xyz_z_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_z_yy[i] = 4.0 * g_xyz_z_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_z_yz[i] = 4.0 * g_xyz_z_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_z_zz[i] = 4.0 * g_xyz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xxz_x_x_xx, g_xxz_x_x_xy, g_xxz_x_x_xz, g_xxz_x_x_yy, g_xxz_x_x_yz, g_xxz_x_x_zz, g_xz_0_0_0_x_x_x_xx, g_xz_0_0_0_x_x_x_xy, g_xz_0_0_0_x_x_x_xz, g_xz_0_0_0_x_x_x_yy, g_xz_0_0_0_x_x_x_yz, g_xz_0_0_0_x_x_x_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_x_xx[i] = -2.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_xxz_x_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_x_xy[i] = -2.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_xxz_x_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_x_xz[i] = -2.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_xxz_x_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_x_yy[i] = -2.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_xxz_x_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_x_yz[i] = -2.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_xxz_x_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_x_zz[i] = -2.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_xxz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xxz_x_y_xx, g_xxz_x_y_xy, g_xxz_x_y_xz, g_xxz_x_y_yy, g_xxz_x_y_yz, g_xxz_x_y_zz, g_xz_0_0_0_x_x_y_xx, g_xz_0_0_0_x_x_y_xy, g_xz_0_0_0_x_x_y_xz, g_xz_0_0_0_x_x_y_yy, g_xz_0_0_0_x_x_y_yz, g_xz_0_0_0_x_x_y_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_y_xx[i] = -2.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_xxz_x_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_y_xy[i] = -2.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_xxz_x_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_y_xz[i] = -2.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_xxz_x_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_y_yy[i] = -2.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_xxz_x_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_y_yz[i] = -2.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_xxz_x_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_y_zz[i] = -2.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_xxz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xxz_x_z_xx, g_xxz_x_z_xy, g_xxz_x_z_xz, g_xxz_x_z_yy, g_xxz_x_z_yz, g_xxz_x_z_zz, g_xz_0_0_0_x_x_z_xx, g_xz_0_0_0_x_x_z_xy, g_xz_0_0_0_x_x_z_xz, g_xz_0_0_0_x_x_z_yy, g_xz_0_0_0_x_x_z_yz, g_xz_0_0_0_x_x_z_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_z_xx[i] = -2.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_xxz_x_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_z_xy[i] = -2.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_xxz_x_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_z_xz[i] = -2.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_xxz_x_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_z_yy[i] = -2.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_xxz_x_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_z_yz[i] = -2.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_xxz_x_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_z_zz[i] = -2.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_xxz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_xxz_y_x_xx, g_xxz_y_x_xy, g_xxz_y_x_xz, g_xxz_y_x_yy, g_xxz_y_x_yz, g_xxz_y_x_zz, g_xz_0_0_0_x_y_x_xx, g_xz_0_0_0_x_y_x_xy, g_xz_0_0_0_x_y_x_xz, g_xz_0_0_0_x_y_x_yy, g_xz_0_0_0_x_y_x_yz, g_xz_0_0_0_x_y_x_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_x_xx[i] = -2.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_xxz_y_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_x_xy[i] = -2.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_xxz_y_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_x_xz[i] = -2.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_xxz_y_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_x_yy[i] = -2.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_xxz_y_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_x_yz[i] = -2.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_xxz_y_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_x_zz[i] = -2.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_xxz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_xxz_y_y_xx, g_xxz_y_y_xy, g_xxz_y_y_xz, g_xxz_y_y_yy, g_xxz_y_y_yz, g_xxz_y_y_zz, g_xz_0_0_0_x_y_y_xx, g_xz_0_0_0_x_y_y_xy, g_xz_0_0_0_x_y_y_xz, g_xz_0_0_0_x_y_y_yy, g_xz_0_0_0_x_y_y_yz, g_xz_0_0_0_x_y_y_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_y_xx[i] = -2.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_xxz_y_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_y_xy[i] = -2.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_xxz_y_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_y_xz[i] = -2.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_xxz_y_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_y_yy[i] = -2.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_xxz_y_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_y_yz[i] = -2.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_xxz_y_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_y_zz[i] = -2.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_xxz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_xxz_y_z_xx, g_xxz_y_z_xy, g_xxz_y_z_xz, g_xxz_y_z_yy, g_xxz_y_z_yz, g_xxz_y_z_zz, g_xz_0_0_0_x_y_z_xx, g_xz_0_0_0_x_y_z_xy, g_xz_0_0_0_x_y_z_xz, g_xz_0_0_0_x_y_z_yy, g_xz_0_0_0_x_y_z_yz, g_xz_0_0_0_x_y_z_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_z_xx[i] = -2.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_xxz_y_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_z_xy[i] = -2.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_xxz_y_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_z_xz[i] = -2.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_xxz_y_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_z_yy[i] = -2.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_xxz_y_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_z_yz[i] = -2.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_xxz_y_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_z_zz[i] = -2.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_xxz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_xxz_z_x_xx, g_xxz_z_x_xy, g_xxz_z_x_xz, g_xxz_z_x_yy, g_xxz_z_x_yz, g_xxz_z_x_zz, g_xz_0_0_0_x_z_x_xx, g_xz_0_0_0_x_z_x_xy, g_xz_0_0_0_x_z_x_xz, g_xz_0_0_0_x_z_x_yy, g_xz_0_0_0_x_z_x_yz, g_xz_0_0_0_x_z_x_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_x_xx[i] = -2.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_xxz_z_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_x_xy[i] = -2.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_xxz_z_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_x_xz[i] = -2.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_xxz_z_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_x_yy[i] = -2.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_xxz_z_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_x_yz[i] = -2.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_xxz_z_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_x_zz[i] = -2.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_xxz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_xxz_z_y_xx, g_xxz_z_y_xy, g_xxz_z_y_xz, g_xxz_z_y_yy, g_xxz_z_y_yz, g_xxz_z_y_zz, g_xz_0_0_0_x_z_y_xx, g_xz_0_0_0_x_z_y_xy, g_xz_0_0_0_x_z_y_xz, g_xz_0_0_0_x_z_y_yy, g_xz_0_0_0_x_z_y_yz, g_xz_0_0_0_x_z_y_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_y_xx[i] = -2.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_xxz_z_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_y_xy[i] = -2.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_xxz_z_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_y_xz[i] = -2.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_xxz_z_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_y_yy[i] = -2.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_xxz_z_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_y_yz[i] = -2.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_xxz_z_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_y_zz[i] = -2.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_xxz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_xxz_z_z_xx, g_xxz_z_z_xy, g_xxz_z_z_xz, g_xxz_z_z_yy, g_xxz_z_z_yz, g_xxz_z_z_zz, g_xz_0_0_0_x_z_z_xx, g_xz_0_0_0_x_z_z_xy, g_xz_0_0_0_x_z_z_xz, g_xz_0_0_0_x_z_z_yy, g_xz_0_0_0_x_z_z_yz, g_xz_0_0_0_x_z_z_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_z_xx[i] = -2.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_xxz_z_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_z_xy[i] = -2.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_xxz_z_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_z_xz[i] = -2.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_xxz_z_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_z_yy[i] = -2.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_xxz_z_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_z_yz[i] = -2.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_xxz_z_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_z_zz[i] = -2.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_xxz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_xyz_x_x_xx, g_xyz_x_x_xy, g_xyz_x_x_xz, g_xyz_x_x_yy, g_xyz_x_x_yz, g_xyz_x_x_zz, g_xz_0_0_0_y_x_x_xx, g_xz_0_0_0_y_x_x_xy, g_xz_0_0_0_y_x_x_xz, g_xz_0_0_0_y_x_x_yy, g_xz_0_0_0_y_x_x_yz, g_xz_0_0_0_y_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_x_xx[i] = 4.0 * g_xyz_x_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_x_xy[i] = 4.0 * g_xyz_x_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_x_xz[i] = 4.0 * g_xyz_x_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_x_yy[i] = 4.0 * g_xyz_x_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_x_yz[i] = 4.0 * g_xyz_x_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_x_zz[i] = 4.0 * g_xyz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_xyz_x_y_xx, g_xyz_x_y_xy, g_xyz_x_y_xz, g_xyz_x_y_yy, g_xyz_x_y_yz, g_xyz_x_y_zz, g_xz_0_0_0_y_x_y_xx, g_xz_0_0_0_y_x_y_xy, g_xz_0_0_0_y_x_y_xz, g_xz_0_0_0_y_x_y_yy, g_xz_0_0_0_y_x_y_yz, g_xz_0_0_0_y_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_y_xx[i] = 4.0 * g_xyz_x_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_y_xy[i] = 4.0 * g_xyz_x_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_y_xz[i] = 4.0 * g_xyz_x_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_y_yy[i] = 4.0 * g_xyz_x_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_y_yz[i] = 4.0 * g_xyz_x_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_y_zz[i] = 4.0 * g_xyz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_xyz_x_z_xx, g_xyz_x_z_xy, g_xyz_x_z_xz, g_xyz_x_z_yy, g_xyz_x_z_yz, g_xyz_x_z_zz, g_xz_0_0_0_y_x_z_xx, g_xz_0_0_0_y_x_z_xy, g_xz_0_0_0_y_x_z_xz, g_xz_0_0_0_y_x_z_yy, g_xz_0_0_0_y_x_z_yz, g_xz_0_0_0_y_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_z_xx[i] = 4.0 * g_xyz_x_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_z_xy[i] = 4.0 * g_xyz_x_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_z_xz[i] = 4.0 * g_xyz_x_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_z_yy[i] = 4.0 * g_xyz_x_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_z_yz[i] = 4.0 * g_xyz_x_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_z_zz[i] = 4.0 * g_xyz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_xyz_y_x_xx, g_xyz_y_x_xy, g_xyz_y_x_xz, g_xyz_y_x_yy, g_xyz_y_x_yz, g_xyz_y_x_zz, g_xz_0_0_0_y_y_x_xx, g_xz_0_0_0_y_y_x_xy, g_xz_0_0_0_y_y_x_xz, g_xz_0_0_0_y_y_x_yy, g_xz_0_0_0_y_y_x_yz, g_xz_0_0_0_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_x_xx[i] = 4.0 * g_xyz_y_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_x_xy[i] = 4.0 * g_xyz_y_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_x_xz[i] = 4.0 * g_xyz_y_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_x_yy[i] = 4.0 * g_xyz_y_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_x_yz[i] = 4.0 * g_xyz_y_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_x_zz[i] = 4.0 * g_xyz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_xyz_y_y_xx, g_xyz_y_y_xy, g_xyz_y_y_xz, g_xyz_y_y_yy, g_xyz_y_y_yz, g_xyz_y_y_zz, g_xz_0_0_0_y_y_y_xx, g_xz_0_0_0_y_y_y_xy, g_xz_0_0_0_y_y_y_xz, g_xz_0_0_0_y_y_y_yy, g_xz_0_0_0_y_y_y_yz, g_xz_0_0_0_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_y_xx[i] = 4.0 * g_xyz_y_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_y_xy[i] = 4.0 * g_xyz_y_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_y_xz[i] = 4.0 * g_xyz_y_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_y_yy[i] = 4.0 * g_xyz_y_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_y_yz[i] = 4.0 * g_xyz_y_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_y_zz[i] = 4.0 * g_xyz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_xyz_y_z_xx, g_xyz_y_z_xy, g_xyz_y_z_xz, g_xyz_y_z_yy, g_xyz_y_z_yz, g_xyz_y_z_zz, g_xz_0_0_0_y_y_z_xx, g_xz_0_0_0_y_y_z_xy, g_xz_0_0_0_y_y_z_xz, g_xz_0_0_0_y_y_z_yy, g_xz_0_0_0_y_y_z_yz, g_xz_0_0_0_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_z_xx[i] = 4.0 * g_xyz_y_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_z_xy[i] = 4.0 * g_xyz_y_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_z_xz[i] = 4.0 * g_xyz_y_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_z_yy[i] = 4.0 * g_xyz_y_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_z_yz[i] = 4.0 * g_xyz_y_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_z_zz[i] = 4.0 * g_xyz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_xyz_z_x_xx, g_xyz_z_x_xy, g_xyz_z_x_xz, g_xyz_z_x_yy, g_xyz_z_x_yz, g_xyz_z_x_zz, g_xz_0_0_0_y_z_x_xx, g_xz_0_0_0_y_z_x_xy, g_xz_0_0_0_y_z_x_xz, g_xz_0_0_0_y_z_x_yy, g_xz_0_0_0_y_z_x_yz, g_xz_0_0_0_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_x_xx[i] = 4.0 * g_xyz_z_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_x_xy[i] = 4.0 * g_xyz_z_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_x_xz[i] = 4.0 * g_xyz_z_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_x_yy[i] = 4.0 * g_xyz_z_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_x_yz[i] = 4.0 * g_xyz_z_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_x_zz[i] = 4.0 * g_xyz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_xyz_z_y_xx, g_xyz_z_y_xy, g_xyz_z_y_xz, g_xyz_z_y_yy, g_xyz_z_y_yz, g_xyz_z_y_zz, g_xz_0_0_0_y_z_y_xx, g_xz_0_0_0_y_z_y_xy, g_xz_0_0_0_y_z_y_xz, g_xz_0_0_0_y_z_y_yy, g_xz_0_0_0_y_z_y_yz, g_xz_0_0_0_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_y_xx[i] = 4.0 * g_xyz_z_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_y_xy[i] = 4.0 * g_xyz_z_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_y_xz[i] = 4.0 * g_xyz_z_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_y_yy[i] = 4.0 * g_xyz_z_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_y_yz[i] = 4.0 * g_xyz_z_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_y_zz[i] = 4.0 * g_xyz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_xyz_z_z_xx, g_xyz_z_z_xy, g_xyz_z_z_xz, g_xyz_z_z_yy, g_xyz_z_z_yz, g_xyz_z_z_zz, g_xz_0_0_0_y_z_z_xx, g_xz_0_0_0_y_z_z_xy, g_xz_0_0_0_y_z_z_xz, g_xz_0_0_0_y_z_z_yy, g_xz_0_0_0_y_z_z_yz, g_xz_0_0_0_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_z_xx[i] = 4.0 * g_xyz_z_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_z_xy[i] = 4.0 * g_xyz_z_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_z_xz[i] = 4.0 * g_xyz_z_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_z_yy[i] = 4.0 * g_xyz_z_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_z_yz[i] = 4.0 * g_xyz_z_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_z_zz[i] = 4.0 * g_xyz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_xz_0_0_0_z_x_x_xx, g_xz_0_0_0_z_x_x_xy, g_xz_0_0_0_z_x_x_xz, g_xz_0_0_0_z_x_x_yy, g_xz_0_0_0_z_x_x_yz, g_xz_0_0_0_z_x_x_zz, g_xzz_x_x_xx, g_xzz_x_x_xy, g_xzz_x_x_xz, g_xzz_x_x_yy, g_xzz_x_x_yz, g_xzz_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_x_xx[i] = -2.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_xzz_x_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_x_xy[i] = -2.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_xzz_x_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_x_xz[i] = -2.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_xzz_x_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_x_yy[i] = -2.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_xzz_x_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_x_yz[i] = -2.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_xzz_x_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_x_zz[i] = -2.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_xzz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, g_xz_0_0_0_z_x_y_xx, g_xz_0_0_0_z_x_y_xy, g_xz_0_0_0_z_x_y_xz, g_xz_0_0_0_z_x_y_yy, g_xz_0_0_0_z_x_y_yz, g_xz_0_0_0_z_x_y_zz, g_xzz_x_y_xx, g_xzz_x_y_xy, g_xzz_x_y_xz, g_xzz_x_y_yy, g_xzz_x_y_yz, g_xzz_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_y_xx[i] = -2.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_xzz_x_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_y_xy[i] = -2.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_xzz_x_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_y_xz[i] = -2.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_xzz_x_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_y_yy[i] = -2.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_xzz_x_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_y_yz[i] = -2.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_xzz_x_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_y_zz[i] = -2.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_xzz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, g_xz_0_0_0_z_x_z_xx, g_xz_0_0_0_z_x_z_xy, g_xz_0_0_0_z_x_z_xz, g_xz_0_0_0_z_x_z_yy, g_xz_0_0_0_z_x_z_yz, g_xz_0_0_0_z_x_z_zz, g_xzz_x_z_xx, g_xzz_x_z_xy, g_xzz_x_z_xz, g_xzz_x_z_yy, g_xzz_x_z_yz, g_xzz_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_z_xx[i] = -2.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_xzz_x_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_z_xy[i] = -2.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_xzz_x_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_z_xz[i] = -2.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_xzz_x_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_z_yy[i] = -2.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_xzz_x_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_z_yz[i] = -2.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_xzz_x_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_z_zz[i] = -2.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_xzz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_xz_0_0_0_z_y_x_xx, g_xz_0_0_0_z_y_x_xy, g_xz_0_0_0_z_y_x_xz, g_xz_0_0_0_z_y_x_yy, g_xz_0_0_0_z_y_x_yz, g_xz_0_0_0_z_y_x_zz, g_xzz_y_x_xx, g_xzz_y_x_xy, g_xzz_y_x_xz, g_xzz_y_x_yy, g_xzz_y_x_yz, g_xzz_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_x_xx[i] = -2.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_xzz_y_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_x_xy[i] = -2.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_xzz_y_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_x_xz[i] = -2.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_xzz_y_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_x_yy[i] = -2.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_xzz_y_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_x_yz[i] = -2.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_xzz_y_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_x_zz[i] = -2.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_xzz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, g_xz_0_0_0_z_y_y_xx, g_xz_0_0_0_z_y_y_xy, g_xz_0_0_0_z_y_y_xz, g_xz_0_0_0_z_y_y_yy, g_xz_0_0_0_z_y_y_yz, g_xz_0_0_0_z_y_y_zz, g_xzz_y_y_xx, g_xzz_y_y_xy, g_xzz_y_y_xz, g_xzz_y_y_yy, g_xzz_y_y_yz, g_xzz_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_y_xx[i] = -2.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_xzz_y_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_y_xy[i] = -2.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_xzz_y_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_y_xz[i] = -2.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_xzz_y_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_y_yy[i] = -2.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_xzz_y_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_y_yz[i] = -2.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_xzz_y_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_y_zz[i] = -2.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_xzz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, g_xz_0_0_0_z_y_z_xx, g_xz_0_0_0_z_y_z_xy, g_xz_0_0_0_z_y_z_xz, g_xz_0_0_0_z_y_z_yy, g_xz_0_0_0_z_y_z_yz, g_xz_0_0_0_z_y_z_zz, g_xzz_y_z_xx, g_xzz_y_z_xy, g_xzz_y_z_xz, g_xzz_y_z_yy, g_xzz_y_z_yz, g_xzz_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_z_xx[i] = -2.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_xzz_y_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_z_xy[i] = -2.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_xzz_y_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_z_xz[i] = -2.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_xzz_y_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_z_yy[i] = -2.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_xzz_y_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_z_yz[i] = -2.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_xzz_y_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_z_zz[i] = -2.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_xzz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, g_xz_0_0_0_z_z_x_xx, g_xz_0_0_0_z_z_x_xy, g_xz_0_0_0_z_z_x_xz, g_xz_0_0_0_z_z_x_yy, g_xz_0_0_0_z_z_x_yz, g_xz_0_0_0_z_z_x_zz, g_xzz_z_x_xx, g_xzz_z_x_xy, g_xzz_z_x_xz, g_xzz_z_x_yy, g_xzz_z_x_yz, g_xzz_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_x_xx[i] = -2.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_xzz_z_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_x_xy[i] = -2.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_xzz_z_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_x_xz[i] = -2.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_xzz_z_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_x_yy[i] = -2.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_xzz_z_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_x_yz[i] = -2.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_xzz_z_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_x_zz[i] = -2.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_xzz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz, g_xz_0_0_0_z_z_y_xx, g_xz_0_0_0_z_z_y_xy, g_xz_0_0_0_z_z_y_xz, g_xz_0_0_0_z_z_y_yy, g_xz_0_0_0_z_z_y_yz, g_xz_0_0_0_z_z_y_zz, g_xzz_z_y_xx, g_xzz_z_y_xy, g_xzz_z_y_xz, g_xzz_z_y_yy, g_xzz_z_y_yz, g_xzz_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_y_xx[i] = -2.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_xzz_z_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_y_xy[i] = -2.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_xzz_z_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_y_xz[i] = -2.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_xzz_z_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_y_yy[i] = -2.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_xzz_z_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_y_yz[i] = -2.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_xzz_z_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_y_zz[i] = -2.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_xzz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz, g_xz_0_0_0_z_z_z_xx, g_xz_0_0_0_z_z_z_xy, g_xz_0_0_0_z_z_z_xz, g_xz_0_0_0_z_z_z_yy, g_xz_0_0_0_z_z_z_yz, g_xz_0_0_0_z_z_z_zz, g_xzz_z_z_xx, g_xzz_z_z_xy, g_xzz_z_z_xz, g_xzz_z_z_yy, g_xzz_z_z_yz, g_xzz_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_z_xx[i] = -2.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_xzz_z_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_z_xy[i] = -2.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_xzz_z_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_z_xz[i] = -2.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_xzz_z_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_z_yy[i] = -2.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_xzz_z_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_z_yz[i] = -2.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_xzz_z_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_z_zz[i] = -2.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_xzz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_xyy_x_x_xx, g_xyy_x_x_xy, g_xyy_x_x_xz, g_xyy_x_x_yy, g_xyy_x_x_yz, g_xyy_x_x_zz, g_yy_0_0_0_x_x_x_xx, g_yy_0_0_0_x_x_x_xy, g_yy_0_0_0_x_x_x_xz, g_yy_0_0_0_x_x_x_yy, g_yy_0_0_0_x_x_x_yz, g_yy_0_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_x_xx[i] = -2.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_xyy_x_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_x_xy[i] = -2.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_xyy_x_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_x_xz[i] = -2.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_xyy_x_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_x_yy[i] = -2.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_xyy_x_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_x_yz[i] = -2.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_xyy_x_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_x_zz[i] = -2.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_xyy_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, g_xyy_x_y_xx, g_xyy_x_y_xy, g_xyy_x_y_xz, g_xyy_x_y_yy, g_xyy_x_y_yz, g_xyy_x_y_zz, g_yy_0_0_0_x_x_y_xx, g_yy_0_0_0_x_x_y_xy, g_yy_0_0_0_x_x_y_xz, g_yy_0_0_0_x_x_y_yy, g_yy_0_0_0_x_x_y_yz, g_yy_0_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_y_xx[i] = -2.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_xyy_x_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_y_xy[i] = -2.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_xyy_x_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_y_xz[i] = -2.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_xyy_x_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_y_yy[i] = -2.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_xyy_x_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_y_yz[i] = -2.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_xyy_x_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_y_zz[i] = -2.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_xyy_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, g_xyy_x_z_xx, g_xyy_x_z_xy, g_xyy_x_z_xz, g_xyy_x_z_yy, g_xyy_x_z_yz, g_xyy_x_z_zz, g_yy_0_0_0_x_x_z_xx, g_yy_0_0_0_x_x_z_xy, g_yy_0_0_0_x_x_z_xz, g_yy_0_0_0_x_x_z_yy, g_yy_0_0_0_x_x_z_yz, g_yy_0_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_z_xx[i] = -2.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_xyy_x_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_z_xy[i] = -2.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_xyy_x_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_z_xz[i] = -2.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_xyy_x_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_z_yy[i] = -2.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_xyy_x_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_z_yz[i] = -2.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_xyy_x_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_z_zz[i] = -2.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_xyy_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_xyy_y_x_xx, g_xyy_y_x_xy, g_xyy_y_x_xz, g_xyy_y_x_yy, g_xyy_y_x_yz, g_xyy_y_x_zz, g_yy_0_0_0_x_y_x_xx, g_yy_0_0_0_x_y_x_xy, g_yy_0_0_0_x_y_x_xz, g_yy_0_0_0_x_y_x_yy, g_yy_0_0_0_x_y_x_yz, g_yy_0_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_x_xx[i] = -2.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_xyy_y_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_x_xy[i] = -2.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_xyy_y_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_x_xz[i] = -2.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_xyy_y_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_x_yy[i] = -2.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_xyy_y_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_x_yz[i] = -2.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_xyy_y_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_x_zz[i] = -2.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_xyy_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, g_xyy_y_y_xx, g_xyy_y_y_xy, g_xyy_y_y_xz, g_xyy_y_y_yy, g_xyy_y_y_yz, g_xyy_y_y_zz, g_yy_0_0_0_x_y_y_xx, g_yy_0_0_0_x_y_y_xy, g_yy_0_0_0_x_y_y_xz, g_yy_0_0_0_x_y_y_yy, g_yy_0_0_0_x_y_y_yz, g_yy_0_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_y_xx[i] = -2.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_xyy_y_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_y_xy[i] = -2.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_xyy_y_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_y_xz[i] = -2.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_xyy_y_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_y_yy[i] = -2.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_xyy_y_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_y_yz[i] = -2.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_xyy_y_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_y_zz[i] = -2.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_xyy_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, g_xyy_y_z_xx, g_xyy_y_z_xy, g_xyy_y_z_xz, g_xyy_y_z_yy, g_xyy_y_z_yz, g_xyy_y_z_zz, g_yy_0_0_0_x_y_z_xx, g_yy_0_0_0_x_y_z_xy, g_yy_0_0_0_x_y_z_xz, g_yy_0_0_0_x_y_z_yy, g_yy_0_0_0_x_y_z_yz, g_yy_0_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_z_xx[i] = -2.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_xyy_y_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_z_xy[i] = -2.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_xyy_y_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_z_xz[i] = -2.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_xyy_y_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_z_yy[i] = -2.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_xyy_y_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_z_yz[i] = -2.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_xyy_y_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_z_zz[i] = -2.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_xyy_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, g_xyy_z_x_xx, g_xyy_z_x_xy, g_xyy_z_x_xz, g_xyy_z_x_yy, g_xyy_z_x_yz, g_xyy_z_x_zz, g_yy_0_0_0_x_z_x_xx, g_yy_0_0_0_x_z_x_xy, g_yy_0_0_0_x_z_x_xz, g_yy_0_0_0_x_z_x_yy, g_yy_0_0_0_x_z_x_yz, g_yy_0_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_x_xx[i] = -2.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_xyy_z_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_x_xy[i] = -2.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_xyy_z_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_x_xz[i] = -2.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_xyy_z_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_x_yy[i] = -2.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_xyy_z_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_x_yz[i] = -2.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_xyy_z_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_x_zz[i] = -2.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_xyy_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz, g_xyy_z_y_xx, g_xyy_z_y_xy, g_xyy_z_y_xz, g_xyy_z_y_yy, g_xyy_z_y_yz, g_xyy_z_y_zz, g_yy_0_0_0_x_z_y_xx, g_yy_0_0_0_x_z_y_xy, g_yy_0_0_0_x_z_y_xz, g_yy_0_0_0_x_z_y_yy, g_yy_0_0_0_x_z_y_yz, g_yy_0_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_y_xx[i] = -2.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_xyy_z_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_y_xy[i] = -2.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_xyy_z_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_y_xz[i] = -2.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_xyy_z_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_y_yy[i] = -2.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_xyy_z_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_y_yz[i] = -2.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_xyy_z_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_y_zz[i] = -2.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_xyy_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz, g_xyy_z_z_xx, g_xyy_z_z_xy, g_xyy_z_z_xz, g_xyy_z_z_yy, g_xyy_z_z_yz, g_xyy_z_z_zz, g_yy_0_0_0_x_z_z_xx, g_yy_0_0_0_x_z_z_xy, g_yy_0_0_0_x_z_z_xz, g_yy_0_0_0_x_z_z_yy, g_yy_0_0_0_x_z_z_yz, g_yy_0_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_z_xx[i] = -2.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_xyy_z_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_z_xy[i] = -2.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_xyy_z_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_z_xz[i] = -2.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_xyy_z_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_z_yy[i] = -2.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_xyy_z_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_z_yz[i] = -2.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_xyy_z_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_z_zz[i] = -2.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_xyy_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, g_yy_0_0_0_y_x_x_xx, g_yy_0_0_0_y_x_x_xy, g_yy_0_0_0_y_x_x_xz, g_yy_0_0_0_y_x_x_yy, g_yy_0_0_0_y_x_x_yz, g_yy_0_0_0_y_x_x_zz, g_yyy_x_x_xx, g_yyy_x_x_xy, g_yyy_x_x_xz, g_yyy_x_x_yy, g_yyy_x_x_yz, g_yyy_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_x_xx[i] = -6.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_yyy_x_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_x_xy[i] = -6.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_yyy_x_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_x_xz[i] = -6.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_yyy_x_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_x_yy[i] = -6.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_yyy_x_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_x_yz[i] = -6.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_yyy_x_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_x_zz[i] = -6.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_yyy_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz, g_yy_0_0_0_y_x_y_xx, g_yy_0_0_0_y_x_y_xy, g_yy_0_0_0_y_x_y_xz, g_yy_0_0_0_y_x_y_yy, g_yy_0_0_0_y_x_y_yz, g_yy_0_0_0_y_x_y_zz, g_yyy_x_y_xx, g_yyy_x_y_xy, g_yyy_x_y_xz, g_yyy_x_y_yy, g_yyy_x_y_yz, g_yyy_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_y_xx[i] = -6.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_yyy_x_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_y_xy[i] = -6.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_yyy_x_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_y_xz[i] = -6.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_yyy_x_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_y_yy[i] = -6.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_yyy_x_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_y_yz[i] = -6.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_yyy_x_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_y_zz[i] = -6.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_yyy_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz, g_yy_0_0_0_y_x_z_xx, g_yy_0_0_0_y_x_z_xy, g_yy_0_0_0_y_x_z_xz, g_yy_0_0_0_y_x_z_yy, g_yy_0_0_0_y_x_z_yz, g_yy_0_0_0_y_x_z_zz, g_yyy_x_z_xx, g_yyy_x_z_xy, g_yyy_x_z_xz, g_yyy_x_z_yy, g_yyy_x_z_yz, g_yyy_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_z_xx[i] = -6.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_yyy_x_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_z_xy[i] = -6.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_yyy_x_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_z_xz[i] = -6.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_yyy_x_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_z_yy[i] = -6.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_yyy_x_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_z_yz[i] = -6.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_yyy_x_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_z_zz[i] = -6.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_yyy_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz, g_yy_0_0_0_y_y_x_xx, g_yy_0_0_0_y_y_x_xy, g_yy_0_0_0_y_y_x_xz, g_yy_0_0_0_y_y_x_yy, g_yy_0_0_0_y_y_x_yz, g_yy_0_0_0_y_y_x_zz, g_yyy_y_x_xx, g_yyy_y_x_xy, g_yyy_y_x_xz, g_yyy_y_x_yy, g_yyy_y_x_yz, g_yyy_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_x_xx[i] = -6.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_yyy_y_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_x_xy[i] = -6.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_yyy_y_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_x_xz[i] = -6.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_yyy_y_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_x_yy[i] = -6.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_yyy_y_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_x_yz[i] = -6.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_yyy_y_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_x_zz[i] = -6.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_yyy_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz, g_yy_0_0_0_y_y_y_xx, g_yy_0_0_0_y_y_y_xy, g_yy_0_0_0_y_y_y_xz, g_yy_0_0_0_y_y_y_yy, g_yy_0_0_0_y_y_y_yz, g_yy_0_0_0_y_y_y_zz, g_yyy_y_y_xx, g_yyy_y_y_xy, g_yyy_y_y_xz, g_yyy_y_y_yy, g_yyy_y_y_yz, g_yyy_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_y_xx[i] = -6.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_yyy_y_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_y_xy[i] = -6.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_yyy_y_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_y_xz[i] = -6.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_yyy_y_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_y_yy[i] = -6.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_yyy_y_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_y_yz[i] = -6.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_yyy_y_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_y_zz[i] = -6.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_yyy_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz, g_yy_0_0_0_y_y_z_xx, g_yy_0_0_0_y_y_z_xy, g_yy_0_0_0_y_y_z_xz, g_yy_0_0_0_y_y_z_yy, g_yy_0_0_0_y_y_z_yz, g_yy_0_0_0_y_y_z_zz, g_yyy_y_z_xx, g_yyy_y_z_xy, g_yyy_y_z_xz, g_yyy_y_z_yy, g_yyy_y_z_yz, g_yyy_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_z_xx[i] = -6.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_yyy_y_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_z_xy[i] = -6.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_yyy_y_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_z_xz[i] = -6.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_yyy_y_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_z_yy[i] = -6.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_yyy_y_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_z_yz[i] = -6.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_yyy_y_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_z_zz[i] = -6.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_yyy_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz, g_yy_0_0_0_y_z_x_xx, g_yy_0_0_0_y_z_x_xy, g_yy_0_0_0_y_z_x_xz, g_yy_0_0_0_y_z_x_yy, g_yy_0_0_0_y_z_x_yz, g_yy_0_0_0_y_z_x_zz, g_yyy_z_x_xx, g_yyy_z_x_xy, g_yyy_z_x_xz, g_yyy_z_x_yy, g_yyy_z_x_yz, g_yyy_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_x_xx[i] = -6.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_yyy_z_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_x_xy[i] = -6.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_yyy_z_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_x_xz[i] = -6.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_yyy_z_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_x_yy[i] = -6.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_yyy_z_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_x_yz[i] = -6.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_yyy_z_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_x_zz[i] = -6.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_yyy_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz, g_yy_0_0_0_y_z_y_xx, g_yy_0_0_0_y_z_y_xy, g_yy_0_0_0_y_z_y_xz, g_yy_0_0_0_y_z_y_yy, g_yy_0_0_0_y_z_y_yz, g_yy_0_0_0_y_z_y_zz, g_yyy_z_y_xx, g_yyy_z_y_xy, g_yyy_z_y_xz, g_yyy_z_y_yy, g_yyy_z_y_yz, g_yyy_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_y_xx[i] = -6.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_yyy_z_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_y_xy[i] = -6.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_yyy_z_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_y_xz[i] = -6.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_yyy_z_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_y_yy[i] = -6.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_yyy_z_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_y_yz[i] = -6.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_yyy_z_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_y_zz[i] = -6.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_yyy_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz, g_yy_0_0_0_y_z_z_xx, g_yy_0_0_0_y_z_z_xy, g_yy_0_0_0_y_z_z_xz, g_yy_0_0_0_y_z_z_yy, g_yy_0_0_0_y_z_z_yz, g_yy_0_0_0_y_z_z_zz, g_yyy_z_z_xx, g_yyy_z_z_xy, g_yyy_z_z_xz, g_yyy_z_z_yy, g_yyy_z_z_yz, g_yyy_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_z_xx[i] = -6.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_yyy_z_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_z_xy[i] = -6.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_yyy_z_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_z_xz[i] = -6.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_yyy_z_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_z_yy[i] = -6.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_yyy_z_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_z_yz[i] = -6.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_yyy_z_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_z_zz[i] = -6.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_yyy_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_x_xx, g_yy_0_0_0_z_x_x_xy, g_yy_0_0_0_z_x_x_xz, g_yy_0_0_0_z_x_x_yy, g_yy_0_0_0_z_x_x_yz, g_yy_0_0_0_z_x_x_zz, g_yyz_x_x_xx, g_yyz_x_x_xy, g_yyz_x_x_xz, g_yyz_x_x_yy, g_yyz_x_x_yz, g_yyz_x_x_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_x_xx[i] = -2.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_yyz_x_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_x_xy[i] = -2.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_yyz_x_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_x_xz[i] = -2.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_yyz_x_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_x_yy[i] = -2.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_yyz_x_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_x_yz[i] = -2.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_yyz_x_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_x_zz[i] = -2.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_yyz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_y_xx, g_yy_0_0_0_z_x_y_xy, g_yy_0_0_0_z_x_y_xz, g_yy_0_0_0_z_x_y_yy, g_yy_0_0_0_z_x_y_yz, g_yy_0_0_0_z_x_y_zz, g_yyz_x_y_xx, g_yyz_x_y_xy, g_yyz_x_y_xz, g_yyz_x_y_yy, g_yyz_x_y_yz, g_yyz_x_y_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_y_xx[i] = -2.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_yyz_x_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_y_xy[i] = -2.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_yyz_x_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_y_xz[i] = -2.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_yyz_x_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_y_yy[i] = -2.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_yyz_x_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_y_yz[i] = -2.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_yyz_x_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_y_zz[i] = -2.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_yyz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_z_xx, g_yy_0_0_0_z_x_z_xy, g_yy_0_0_0_z_x_z_xz, g_yy_0_0_0_z_x_z_yy, g_yy_0_0_0_z_x_z_yz, g_yy_0_0_0_z_x_z_zz, g_yyz_x_z_xx, g_yyz_x_z_xy, g_yyz_x_z_xz, g_yyz_x_z_yy, g_yyz_x_z_yz, g_yyz_x_z_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_z_xx[i] = -2.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_yyz_x_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_z_xy[i] = -2.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_yyz_x_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_z_xz[i] = -2.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_yyz_x_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_z_yy[i] = -2.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_yyz_x_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_z_yz[i] = -2.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_yyz_x_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_z_zz[i] = -2.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_yyz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_x_xx, g_yy_0_0_0_z_y_x_xy, g_yy_0_0_0_z_y_x_xz, g_yy_0_0_0_z_y_x_yy, g_yy_0_0_0_z_y_x_yz, g_yy_0_0_0_z_y_x_zz, g_yyz_y_x_xx, g_yyz_y_x_xy, g_yyz_y_x_xz, g_yyz_y_x_yy, g_yyz_y_x_yz, g_yyz_y_x_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_x_xx[i] = -2.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_yyz_y_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_x_xy[i] = -2.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_yyz_y_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_x_xz[i] = -2.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_yyz_y_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_x_yy[i] = -2.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_yyz_y_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_x_yz[i] = -2.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_yyz_y_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_x_zz[i] = -2.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_yyz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_y_xx, g_yy_0_0_0_z_y_y_xy, g_yy_0_0_0_z_y_y_xz, g_yy_0_0_0_z_y_y_yy, g_yy_0_0_0_z_y_y_yz, g_yy_0_0_0_z_y_y_zz, g_yyz_y_y_xx, g_yyz_y_y_xy, g_yyz_y_y_xz, g_yyz_y_y_yy, g_yyz_y_y_yz, g_yyz_y_y_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_y_xx[i] = -2.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_yyz_y_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_y_xy[i] = -2.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_yyz_y_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_y_xz[i] = -2.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_yyz_y_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_y_yy[i] = -2.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_yyz_y_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_y_yz[i] = -2.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_yyz_y_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_y_zz[i] = -2.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_yyz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_z_xx, g_yy_0_0_0_z_y_z_xy, g_yy_0_0_0_z_y_z_xz, g_yy_0_0_0_z_y_z_yy, g_yy_0_0_0_z_y_z_yz, g_yy_0_0_0_z_y_z_zz, g_yyz_y_z_xx, g_yyz_y_z_xy, g_yyz_y_z_xz, g_yyz_y_z_yy, g_yyz_y_z_yz, g_yyz_y_z_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_z_xx[i] = -2.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_yyz_y_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_z_xy[i] = -2.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_yyz_y_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_z_xz[i] = -2.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_yyz_y_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_z_yy[i] = -2.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_yyz_y_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_z_yz[i] = -2.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_yyz_y_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_z_zz[i] = -2.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_yyz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_x_xx, g_yy_0_0_0_z_z_x_xy, g_yy_0_0_0_z_z_x_xz, g_yy_0_0_0_z_z_x_yy, g_yy_0_0_0_z_z_x_yz, g_yy_0_0_0_z_z_x_zz, g_yyz_z_x_xx, g_yyz_z_x_xy, g_yyz_z_x_xz, g_yyz_z_x_yy, g_yyz_z_x_yz, g_yyz_z_x_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_x_xx[i] = -2.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_yyz_z_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_x_xy[i] = -2.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_yyz_z_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_x_xz[i] = -2.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_yyz_z_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_x_yy[i] = -2.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_yyz_z_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_x_yz[i] = -2.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_yyz_z_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_x_zz[i] = -2.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_yyz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_y_xx, g_yy_0_0_0_z_z_y_xy, g_yy_0_0_0_z_z_y_xz, g_yy_0_0_0_z_z_y_yy, g_yy_0_0_0_z_z_y_yz, g_yy_0_0_0_z_z_y_zz, g_yyz_z_y_xx, g_yyz_z_y_xy, g_yyz_z_y_xz, g_yyz_z_y_yy, g_yyz_z_y_yz, g_yyz_z_y_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_y_xx[i] = -2.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_yyz_z_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_y_xy[i] = -2.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_yyz_z_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_y_xz[i] = -2.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_yyz_z_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_y_yy[i] = -2.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_yyz_z_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_y_yz[i] = -2.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_yyz_z_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_y_zz[i] = -2.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_yyz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_z_xx, g_yy_0_0_0_z_z_z_xy, g_yy_0_0_0_z_z_z_xz, g_yy_0_0_0_z_z_z_yy, g_yy_0_0_0_z_z_z_yz, g_yy_0_0_0_z_z_z_zz, g_yyz_z_z_xx, g_yyz_z_z_xy, g_yyz_z_z_xz, g_yyz_z_z_yy, g_yyz_z_z_yz, g_yyz_z_z_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_z_xx[i] = -2.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_yyz_z_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_z_xy[i] = -2.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_yyz_z_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_z_xz[i] = -2.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_yyz_z_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_z_yy[i] = -2.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_yyz_z_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_z_yz[i] = -2.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_yyz_z_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_z_zz[i] = -2.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_yyz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xyz_x_x_xx, g_xyz_x_x_xy, g_xyz_x_x_xz, g_xyz_x_x_yy, g_xyz_x_x_yz, g_xyz_x_x_zz, g_yz_0_0_0_x_x_x_xx, g_yz_0_0_0_x_x_x_xy, g_yz_0_0_0_x_x_x_xz, g_yz_0_0_0_x_x_x_yy, g_yz_0_0_0_x_x_x_yz, g_yz_0_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_x_xx[i] = 4.0 * g_xyz_x_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_x_xy[i] = 4.0 * g_xyz_x_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_x_xz[i] = 4.0 * g_xyz_x_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_x_yy[i] = 4.0 * g_xyz_x_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_x_yz[i] = 4.0 * g_xyz_x_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_x_zz[i] = 4.0 * g_xyz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xyz_x_y_xx, g_xyz_x_y_xy, g_xyz_x_y_xz, g_xyz_x_y_yy, g_xyz_x_y_yz, g_xyz_x_y_zz, g_yz_0_0_0_x_x_y_xx, g_yz_0_0_0_x_x_y_xy, g_yz_0_0_0_x_x_y_xz, g_yz_0_0_0_x_x_y_yy, g_yz_0_0_0_x_x_y_yz, g_yz_0_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_y_xx[i] = 4.0 * g_xyz_x_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_y_xy[i] = 4.0 * g_xyz_x_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_y_xz[i] = 4.0 * g_xyz_x_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_y_yy[i] = 4.0 * g_xyz_x_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_y_yz[i] = 4.0 * g_xyz_x_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_y_zz[i] = 4.0 * g_xyz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xyz_x_z_xx, g_xyz_x_z_xy, g_xyz_x_z_xz, g_xyz_x_z_yy, g_xyz_x_z_yz, g_xyz_x_z_zz, g_yz_0_0_0_x_x_z_xx, g_yz_0_0_0_x_x_z_xy, g_yz_0_0_0_x_x_z_xz, g_yz_0_0_0_x_x_z_yy, g_yz_0_0_0_x_x_z_yz, g_yz_0_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_z_xx[i] = 4.0 * g_xyz_x_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_z_xy[i] = 4.0 * g_xyz_x_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_z_xz[i] = 4.0 * g_xyz_x_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_z_yy[i] = 4.0 * g_xyz_x_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_z_yz[i] = 4.0 * g_xyz_x_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_z_zz[i] = 4.0 * g_xyz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xyz_y_x_xx, g_xyz_y_x_xy, g_xyz_y_x_xz, g_xyz_y_x_yy, g_xyz_y_x_yz, g_xyz_y_x_zz, g_yz_0_0_0_x_y_x_xx, g_yz_0_0_0_x_y_x_xy, g_yz_0_0_0_x_y_x_xz, g_yz_0_0_0_x_y_x_yy, g_yz_0_0_0_x_y_x_yz, g_yz_0_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_x_xx[i] = 4.0 * g_xyz_y_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_x_xy[i] = 4.0 * g_xyz_y_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_x_xz[i] = 4.0 * g_xyz_y_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_x_yy[i] = 4.0 * g_xyz_y_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_x_yz[i] = 4.0 * g_xyz_y_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_x_zz[i] = 4.0 * g_xyz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xyz_y_y_xx, g_xyz_y_y_xy, g_xyz_y_y_xz, g_xyz_y_y_yy, g_xyz_y_y_yz, g_xyz_y_y_zz, g_yz_0_0_0_x_y_y_xx, g_yz_0_0_0_x_y_y_xy, g_yz_0_0_0_x_y_y_xz, g_yz_0_0_0_x_y_y_yy, g_yz_0_0_0_x_y_y_yz, g_yz_0_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_y_xx[i] = 4.0 * g_xyz_y_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_y_xy[i] = 4.0 * g_xyz_y_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_y_xz[i] = 4.0 * g_xyz_y_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_y_yy[i] = 4.0 * g_xyz_y_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_y_yz[i] = 4.0 * g_xyz_y_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_y_zz[i] = 4.0 * g_xyz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xyz_y_z_xx, g_xyz_y_z_xy, g_xyz_y_z_xz, g_xyz_y_z_yy, g_xyz_y_z_yz, g_xyz_y_z_zz, g_yz_0_0_0_x_y_z_xx, g_yz_0_0_0_x_y_z_xy, g_yz_0_0_0_x_y_z_xz, g_yz_0_0_0_x_y_z_yy, g_yz_0_0_0_x_y_z_yz, g_yz_0_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_z_xx[i] = 4.0 * g_xyz_y_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_z_xy[i] = 4.0 * g_xyz_y_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_z_xz[i] = 4.0 * g_xyz_y_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_z_yy[i] = 4.0 * g_xyz_y_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_z_yz[i] = 4.0 * g_xyz_y_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_z_zz[i] = 4.0 * g_xyz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_xyz_z_x_xx, g_xyz_z_x_xy, g_xyz_z_x_xz, g_xyz_z_x_yy, g_xyz_z_x_yz, g_xyz_z_x_zz, g_yz_0_0_0_x_z_x_xx, g_yz_0_0_0_x_z_x_xy, g_yz_0_0_0_x_z_x_xz, g_yz_0_0_0_x_z_x_yy, g_yz_0_0_0_x_z_x_yz, g_yz_0_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_x_xx[i] = 4.0 * g_xyz_z_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_x_xy[i] = 4.0 * g_xyz_z_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_x_xz[i] = 4.0 * g_xyz_z_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_x_yy[i] = 4.0 * g_xyz_z_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_x_yz[i] = 4.0 * g_xyz_z_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_x_zz[i] = 4.0 * g_xyz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_xyz_z_y_xx, g_xyz_z_y_xy, g_xyz_z_y_xz, g_xyz_z_y_yy, g_xyz_z_y_yz, g_xyz_z_y_zz, g_yz_0_0_0_x_z_y_xx, g_yz_0_0_0_x_z_y_xy, g_yz_0_0_0_x_z_y_xz, g_yz_0_0_0_x_z_y_yy, g_yz_0_0_0_x_z_y_yz, g_yz_0_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_y_xx[i] = 4.0 * g_xyz_z_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_y_xy[i] = 4.0 * g_xyz_z_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_y_xz[i] = 4.0 * g_xyz_z_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_y_yy[i] = 4.0 * g_xyz_z_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_y_yz[i] = 4.0 * g_xyz_z_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_y_zz[i] = 4.0 * g_xyz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_xyz_z_z_xx, g_xyz_z_z_xy, g_xyz_z_z_xz, g_xyz_z_z_yy, g_xyz_z_z_yz, g_xyz_z_z_zz, g_yz_0_0_0_x_z_z_xx, g_yz_0_0_0_x_z_z_xy, g_yz_0_0_0_x_z_z_xz, g_yz_0_0_0_x_z_z_yy, g_yz_0_0_0_x_z_z_yz, g_yz_0_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_z_xx[i] = 4.0 * g_xyz_z_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_z_xy[i] = 4.0 * g_xyz_z_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_z_xz[i] = 4.0 * g_xyz_z_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_z_yy[i] = 4.0 * g_xyz_z_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_z_yz[i] = 4.0 * g_xyz_z_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_z_zz[i] = 4.0 * g_xyz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_yyz_x_x_xx, g_yyz_x_x_xy, g_yyz_x_x_xz, g_yyz_x_x_yy, g_yyz_x_x_yz, g_yyz_x_x_zz, g_yz_0_0_0_y_x_x_xx, g_yz_0_0_0_y_x_x_xy, g_yz_0_0_0_y_x_x_xz, g_yz_0_0_0_y_x_x_yy, g_yz_0_0_0_y_x_x_yz, g_yz_0_0_0_y_x_x_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_x_xx[i] = -2.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_yyz_x_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_x_xy[i] = -2.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_yyz_x_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_x_xz[i] = -2.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_yyz_x_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_x_yy[i] = -2.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_yyz_x_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_x_yz[i] = -2.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_yyz_x_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_x_zz[i] = -2.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_yyz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_yyz_x_y_xx, g_yyz_x_y_xy, g_yyz_x_y_xz, g_yyz_x_y_yy, g_yyz_x_y_yz, g_yyz_x_y_zz, g_yz_0_0_0_y_x_y_xx, g_yz_0_0_0_y_x_y_xy, g_yz_0_0_0_y_x_y_xz, g_yz_0_0_0_y_x_y_yy, g_yz_0_0_0_y_x_y_yz, g_yz_0_0_0_y_x_y_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_y_xx[i] = -2.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_yyz_x_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_y_xy[i] = -2.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_yyz_x_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_y_xz[i] = -2.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_yyz_x_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_y_yy[i] = -2.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_yyz_x_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_y_yz[i] = -2.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_yyz_x_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_y_zz[i] = -2.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_yyz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_yyz_x_z_xx, g_yyz_x_z_xy, g_yyz_x_z_xz, g_yyz_x_z_yy, g_yyz_x_z_yz, g_yyz_x_z_zz, g_yz_0_0_0_y_x_z_xx, g_yz_0_0_0_y_x_z_xy, g_yz_0_0_0_y_x_z_xz, g_yz_0_0_0_y_x_z_yy, g_yz_0_0_0_y_x_z_yz, g_yz_0_0_0_y_x_z_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_z_xx[i] = -2.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_yyz_x_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_z_xy[i] = -2.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_yyz_x_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_z_xz[i] = -2.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_yyz_x_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_z_yy[i] = -2.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_yyz_x_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_z_yz[i] = -2.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_yyz_x_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_z_zz[i] = -2.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_yyz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_yyz_y_x_xx, g_yyz_y_x_xy, g_yyz_y_x_xz, g_yyz_y_x_yy, g_yyz_y_x_yz, g_yyz_y_x_zz, g_yz_0_0_0_y_y_x_xx, g_yz_0_0_0_y_y_x_xy, g_yz_0_0_0_y_y_x_xz, g_yz_0_0_0_y_y_x_yy, g_yz_0_0_0_y_y_x_yz, g_yz_0_0_0_y_y_x_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_x_xx[i] = -2.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_yyz_y_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_x_xy[i] = -2.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_yyz_y_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_x_xz[i] = -2.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_yyz_y_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_x_yy[i] = -2.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_yyz_y_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_x_yz[i] = -2.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_yyz_y_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_x_zz[i] = -2.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_yyz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_yyz_y_y_xx, g_yyz_y_y_xy, g_yyz_y_y_xz, g_yyz_y_y_yy, g_yyz_y_y_yz, g_yyz_y_y_zz, g_yz_0_0_0_y_y_y_xx, g_yz_0_0_0_y_y_y_xy, g_yz_0_0_0_y_y_y_xz, g_yz_0_0_0_y_y_y_yy, g_yz_0_0_0_y_y_y_yz, g_yz_0_0_0_y_y_y_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_y_xx[i] = -2.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_yyz_y_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_y_xy[i] = -2.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_yyz_y_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_y_xz[i] = -2.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_yyz_y_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_y_yy[i] = -2.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_yyz_y_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_y_yz[i] = -2.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_yyz_y_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_y_zz[i] = -2.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_yyz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_yyz_y_z_xx, g_yyz_y_z_xy, g_yyz_y_z_xz, g_yyz_y_z_yy, g_yyz_y_z_yz, g_yyz_y_z_zz, g_yz_0_0_0_y_y_z_xx, g_yz_0_0_0_y_y_z_xy, g_yz_0_0_0_y_y_z_xz, g_yz_0_0_0_y_y_z_yy, g_yz_0_0_0_y_y_z_yz, g_yz_0_0_0_y_y_z_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_z_xx[i] = -2.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_yyz_y_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_z_xy[i] = -2.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_yyz_y_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_z_xz[i] = -2.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_yyz_y_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_z_yy[i] = -2.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_yyz_y_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_z_yz[i] = -2.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_yyz_y_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_z_zz[i] = -2.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_yyz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_yyz_z_x_xx, g_yyz_z_x_xy, g_yyz_z_x_xz, g_yyz_z_x_yy, g_yyz_z_x_yz, g_yyz_z_x_zz, g_yz_0_0_0_y_z_x_xx, g_yz_0_0_0_y_z_x_xy, g_yz_0_0_0_y_z_x_xz, g_yz_0_0_0_y_z_x_yy, g_yz_0_0_0_y_z_x_yz, g_yz_0_0_0_y_z_x_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_x_xx[i] = -2.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_yyz_z_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_x_xy[i] = -2.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_yyz_z_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_x_xz[i] = -2.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_yyz_z_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_x_yy[i] = -2.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_yyz_z_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_x_yz[i] = -2.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_yyz_z_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_x_zz[i] = -2.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_yyz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_yyz_z_y_xx, g_yyz_z_y_xy, g_yyz_z_y_xz, g_yyz_z_y_yy, g_yyz_z_y_yz, g_yyz_z_y_zz, g_yz_0_0_0_y_z_y_xx, g_yz_0_0_0_y_z_y_xy, g_yz_0_0_0_y_z_y_xz, g_yz_0_0_0_y_z_y_yy, g_yz_0_0_0_y_z_y_yz, g_yz_0_0_0_y_z_y_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_y_xx[i] = -2.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_yyz_z_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_y_xy[i] = -2.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_yyz_z_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_y_xz[i] = -2.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_yyz_z_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_y_yy[i] = -2.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_yyz_z_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_y_yz[i] = -2.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_yyz_z_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_y_zz[i] = -2.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_yyz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_yyz_z_z_xx, g_yyz_z_z_xy, g_yyz_z_z_xz, g_yyz_z_z_yy, g_yyz_z_z_yz, g_yyz_z_z_zz, g_yz_0_0_0_y_z_z_xx, g_yz_0_0_0_y_z_z_xy, g_yz_0_0_0_y_z_z_xz, g_yz_0_0_0_y_z_z_yy, g_yz_0_0_0_y_z_z_yz, g_yz_0_0_0_y_z_z_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_z_xx[i] = -2.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_yyz_z_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_z_xy[i] = -2.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_yyz_z_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_z_xz[i] = -2.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_yyz_z_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_z_yy[i] = -2.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_yyz_z_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_z_yz[i] = -2.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_yyz_z_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_z_zz[i] = -2.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_yyz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, g_yz_0_0_0_z_x_x_xx, g_yz_0_0_0_z_x_x_xy, g_yz_0_0_0_z_x_x_xz, g_yz_0_0_0_z_x_x_yy, g_yz_0_0_0_z_x_x_yz, g_yz_0_0_0_z_x_x_zz, g_yzz_x_x_xx, g_yzz_x_x_xy, g_yzz_x_x_xz, g_yzz_x_x_yy, g_yzz_x_x_yz, g_yzz_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_x_xx[i] = -2.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_yzz_x_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_x_xy[i] = -2.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_yzz_x_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_x_xz[i] = -2.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_yzz_x_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_x_yy[i] = -2.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_yzz_x_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_x_yz[i] = -2.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_yzz_x_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_x_zz[i] = -2.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_yzz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz, g_yz_0_0_0_z_x_y_xx, g_yz_0_0_0_z_x_y_xy, g_yz_0_0_0_z_x_y_xz, g_yz_0_0_0_z_x_y_yy, g_yz_0_0_0_z_x_y_yz, g_yz_0_0_0_z_x_y_zz, g_yzz_x_y_xx, g_yzz_x_y_xy, g_yzz_x_y_xz, g_yzz_x_y_yy, g_yzz_x_y_yz, g_yzz_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_y_xx[i] = -2.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_yzz_x_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_y_xy[i] = -2.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_yzz_x_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_y_xz[i] = -2.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_yzz_x_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_y_yy[i] = -2.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_yzz_x_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_y_yz[i] = -2.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_yzz_x_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_y_zz[i] = -2.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_yzz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz, g_yz_0_0_0_z_x_z_xx, g_yz_0_0_0_z_x_z_xy, g_yz_0_0_0_z_x_z_xz, g_yz_0_0_0_z_x_z_yy, g_yz_0_0_0_z_x_z_yz, g_yz_0_0_0_z_x_z_zz, g_yzz_x_z_xx, g_yzz_x_z_xy, g_yzz_x_z_xz, g_yzz_x_z_yy, g_yzz_x_z_yz, g_yzz_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_z_xx[i] = -2.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_yzz_x_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_z_xy[i] = -2.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_yzz_x_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_z_xz[i] = -2.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_yzz_x_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_z_yy[i] = -2.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_yzz_x_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_z_yz[i] = -2.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_yzz_x_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_z_zz[i] = -2.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_yzz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz, g_yz_0_0_0_z_y_x_xx, g_yz_0_0_0_z_y_x_xy, g_yz_0_0_0_z_y_x_xz, g_yz_0_0_0_z_y_x_yy, g_yz_0_0_0_z_y_x_yz, g_yz_0_0_0_z_y_x_zz, g_yzz_y_x_xx, g_yzz_y_x_xy, g_yzz_y_x_xz, g_yzz_y_x_yy, g_yzz_y_x_yz, g_yzz_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_x_xx[i] = -2.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_yzz_y_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_x_xy[i] = -2.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_yzz_y_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_x_xz[i] = -2.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_yzz_y_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_x_yy[i] = -2.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_yzz_y_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_x_yz[i] = -2.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_yzz_y_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_x_zz[i] = -2.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_yzz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz, g_yz_0_0_0_z_y_y_xx, g_yz_0_0_0_z_y_y_xy, g_yz_0_0_0_z_y_y_xz, g_yz_0_0_0_z_y_y_yy, g_yz_0_0_0_z_y_y_yz, g_yz_0_0_0_z_y_y_zz, g_yzz_y_y_xx, g_yzz_y_y_xy, g_yzz_y_y_xz, g_yzz_y_y_yy, g_yzz_y_y_yz, g_yzz_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_y_xx[i] = -2.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_yzz_y_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_y_xy[i] = -2.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_yzz_y_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_y_xz[i] = -2.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_yzz_y_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_y_yy[i] = -2.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_yzz_y_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_y_yz[i] = -2.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_yzz_y_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_y_zz[i] = -2.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_yzz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz, g_yz_0_0_0_z_y_z_xx, g_yz_0_0_0_z_y_z_xy, g_yz_0_0_0_z_y_z_xz, g_yz_0_0_0_z_y_z_yy, g_yz_0_0_0_z_y_z_yz, g_yz_0_0_0_z_y_z_zz, g_yzz_y_z_xx, g_yzz_y_z_xy, g_yzz_y_z_xz, g_yzz_y_z_yy, g_yzz_y_z_yz, g_yzz_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_z_xx[i] = -2.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_yzz_y_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_z_xy[i] = -2.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_yzz_y_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_z_xz[i] = -2.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_yzz_y_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_z_yy[i] = -2.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_yzz_y_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_z_yz[i] = -2.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_yzz_y_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_z_zz[i] = -2.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_yzz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz, g_yz_0_0_0_z_z_x_xx, g_yz_0_0_0_z_z_x_xy, g_yz_0_0_0_z_z_x_xz, g_yz_0_0_0_z_z_x_yy, g_yz_0_0_0_z_z_x_yz, g_yz_0_0_0_z_z_x_zz, g_yzz_z_x_xx, g_yzz_z_x_xy, g_yzz_z_x_xz, g_yzz_z_x_yy, g_yzz_z_x_yz, g_yzz_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_x_xx[i] = -2.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_yzz_z_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_x_xy[i] = -2.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_yzz_z_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_x_xz[i] = -2.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_yzz_z_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_x_yy[i] = -2.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_yzz_z_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_x_yz[i] = -2.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_yzz_z_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_x_zz[i] = -2.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_yzz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz, g_yz_0_0_0_z_z_y_xx, g_yz_0_0_0_z_z_y_xy, g_yz_0_0_0_z_z_y_xz, g_yz_0_0_0_z_z_y_yy, g_yz_0_0_0_z_z_y_yz, g_yz_0_0_0_z_z_y_zz, g_yzz_z_y_xx, g_yzz_z_y_xy, g_yzz_z_y_xz, g_yzz_z_y_yy, g_yzz_z_y_yz, g_yzz_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_y_xx[i] = -2.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_yzz_z_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_y_xy[i] = -2.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_yzz_z_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_y_xz[i] = -2.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_yzz_z_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_y_yy[i] = -2.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_yzz_z_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_y_yz[i] = -2.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_yzz_z_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_y_zz[i] = -2.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_yzz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz, g_yz_0_0_0_z_z_z_xx, g_yz_0_0_0_z_z_z_xy, g_yz_0_0_0_z_z_z_xz, g_yz_0_0_0_z_z_z_yy, g_yz_0_0_0_z_z_z_yz, g_yz_0_0_0_z_z_z_zz, g_yzz_z_z_xx, g_yzz_z_z_xy, g_yzz_z_z_xz, g_yzz_z_z_yy, g_yzz_z_z_yz, g_yzz_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_z_xx[i] = -2.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_yzz_z_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_z_xy[i] = -2.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_yzz_z_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_z_xz[i] = -2.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_yzz_z_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_z_yy[i] = -2.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_yzz_z_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_z_yz[i] = -2.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_yzz_z_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_z_zz[i] = -2.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_yzz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_xzz_x_x_xx, g_xzz_x_x_xy, g_xzz_x_x_xz, g_xzz_x_x_yy, g_xzz_x_x_yz, g_xzz_x_x_zz, g_zz_0_0_0_x_x_x_xx, g_zz_0_0_0_x_x_x_xy, g_zz_0_0_0_x_x_x_xz, g_zz_0_0_0_x_x_x_yy, g_zz_0_0_0_x_x_x_yz, g_zz_0_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_x_xx[i] = -2.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_xzz_x_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_x_xy[i] = -2.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_xzz_x_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_x_xz[i] = -2.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_xzz_x_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_x_yy[i] = -2.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_xzz_x_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_x_yz[i] = -2.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_xzz_x_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_x_zz[i] = -2.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_xzz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, g_xzz_x_y_xx, g_xzz_x_y_xy, g_xzz_x_y_xz, g_xzz_x_y_yy, g_xzz_x_y_yz, g_xzz_x_y_zz, g_zz_0_0_0_x_x_y_xx, g_zz_0_0_0_x_x_y_xy, g_zz_0_0_0_x_x_y_xz, g_zz_0_0_0_x_x_y_yy, g_zz_0_0_0_x_x_y_yz, g_zz_0_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_y_xx[i] = -2.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_xzz_x_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_y_xy[i] = -2.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_xzz_x_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_y_xz[i] = -2.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_xzz_x_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_y_yy[i] = -2.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_xzz_x_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_y_yz[i] = -2.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_xzz_x_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_y_zz[i] = -2.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_xzz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, g_xzz_x_z_xx, g_xzz_x_z_xy, g_xzz_x_z_xz, g_xzz_x_z_yy, g_xzz_x_z_yz, g_xzz_x_z_zz, g_zz_0_0_0_x_x_z_xx, g_zz_0_0_0_x_x_z_xy, g_zz_0_0_0_x_x_z_xz, g_zz_0_0_0_x_x_z_yy, g_zz_0_0_0_x_x_z_yz, g_zz_0_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_z_xx[i] = -2.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_xzz_x_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_z_xy[i] = -2.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_xzz_x_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_z_xz[i] = -2.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_xzz_x_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_z_yy[i] = -2.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_xzz_x_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_z_yz[i] = -2.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_xzz_x_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_z_zz[i] = -2.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_xzz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_xzz_y_x_xx, g_xzz_y_x_xy, g_xzz_y_x_xz, g_xzz_y_x_yy, g_xzz_y_x_yz, g_xzz_y_x_zz, g_zz_0_0_0_x_y_x_xx, g_zz_0_0_0_x_y_x_xy, g_zz_0_0_0_x_y_x_xz, g_zz_0_0_0_x_y_x_yy, g_zz_0_0_0_x_y_x_yz, g_zz_0_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_x_xx[i] = -2.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_xzz_y_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_x_xy[i] = -2.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_xzz_y_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_x_xz[i] = -2.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_xzz_y_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_x_yy[i] = -2.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_xzz_y_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_x_yz[i] = -2.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_xzz_y_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_x_zz[i] = -2.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_xzz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, g_xzz_y_y_xx, g_xzz_y_y_xy, g_xzz_y_y_xz, g_xzz_y_y_yy, g_xzz_y_y_yz, g_xzz_y_y_zz, g_zz_0_0_0_x_y_y_xx, g_zz_0_0_0_x_y_y_xy, g_zz_0_0_0_x_y_y_xz, g_zz_0_0_0_x_y_y_yy, g_zz_0_0_0_x_y_y_yz, g_zz_0_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_y_xx[i] = -2.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_xzz_y_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_y_xy[i] = -2.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_xzz_y_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_y_xz[i] = -2.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_xzz_y_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_y_yy[i] = -2.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_xzz_y_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_y_yz[i] = -2.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_xzz_y_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_y_zz[i] = -2.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_xzz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, g_xzz_y_z_xx, g_xzz_y_z_xy, g_xzz_y_z_xz, g_xzz_y_z_yy, g_xzz_y_z_yz, g_xzz_y_z_zz, g_zz_0_0_0_x_y_z_xx, g_zz_0_0_0_x_y_z_xy, g_zz_0_0_0_x_y_z_xz, g_zz_0_0_0_x_y_z_yy, g_zz_0_0_0_x_y_z_yz, g_zz_0_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_z_xx[i] = -2.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_xzz_y_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_z_xy[i] = -2.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_xzz_y_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_z_xz[i] = -2.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_xzz_y_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_z_yy[i] = -2.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_xzz_y_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_z_yz[i] = -2.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_xzz_y_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_z_zz[i] = -2.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_xzz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, g_xzz_z_x_xx, g_xzz_z_x_xy, g_xzz_z_x_xz, g_xzz_z_x_yy, g_xzz_z_x_yz, g_xzz_z_x_zz, g_zz_0_0_0_x_z_x_xx, g_zz_0_0_0_x_z_x_xy, g_zz_0_0_0_x_z_x_xz, g_zz_0_0_0_x_z_x_yy, g_zz_0_0_0_x_z_x_yz, g_zz_0_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_x_xx[i] = -2.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_xzz_z_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_x_xy[i] = -2.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_xzz_z_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_x_xz[i] = -2.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_xzz_z_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_x_yy[i] = -2.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_xzz_z_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_x_yz[i] = -2.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_xzz_z_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_x_zz[i] = -2.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_xzz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz, g_xzz_z_y_xx, g_xzz_z_y_xy, g_xzz_z_y_xz, g_xzz_z_y_yy, g_xzz_z_y_yz, g_xzz_z_y_zz, g_zz_0_0_0_x_z_y_xx, g_zz_0_0_0_x_z_y_xy, g_zz_0_0_0_x_z_y_xz, g_zz_0_0_0_x_z_y_yy, g_zz_0_0_0_x_z_y_yz, g_zz_0_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_y_xx[i] = -2.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_xzz_z_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_y_xy[i] = -2.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_xzz_z_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_y_xz[i] = -2.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_xzz_z_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_y_yy[i] = -2.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_xzz_z_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_y_yz[i] = -2.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_xzz_z_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_y_zz[i] = -2.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_xzz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz, g_xzz_z_z_xx, g_xzz_z_z_xy, g_xzz_z_z_xz, g_xzz_z_z_yy, g_xzz_z_z_yz, g_xzz_z_z_zz, g_zz_0_0_0_x_z_z_xx, g_zz_0_0_0_x_z_z_xy, g_zz_0_0_0_x_z_z_xz, g_zz_0_0_0_x_z_z_yy, g_zz_0_0_0_x_z_z_yz, g_zz_0_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_z_xx[i] = -2.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_xzz_z_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_z_xy[i] = -2.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_xzz_z_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_z_xz[i] = -2.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_xzz_z_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_z_yy[i] = -2.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_xzz_z_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_z_yz[i] = -2.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_xzz_z_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_z_zz[i] = -2.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_xzz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, g_yzz_x_x_xx, g_yzz_x_x_xy, g_yzz_x_x_xz, g_yzz_x_x_yy, g_yzz_x_x_yz, g_yzz_x_x_zz, g_zz_0_0_0_y_x_x_xx, g_zz_0_0_0_y_x_x_xy, g_zz_0_0_0_y_x_x_xz, g_zz_0_0_0_y_x_x_yy, g_zz_0_0_0_y_x_x_yz, g_zz_0_0_0_y_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_x_xx[i] = -2.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_yzz_x_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_x_xy[i] = -2.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_yzz_x_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_x_xz[i] = -2.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_yzz_x_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_x_yy[i] = -2.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_yzz_x_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_x_yz[i] = -2.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_yzz_x_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_x_zz[i] = -2.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_yzz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz, g_yzz_x_y_xx, g_yzz_x_y_xy, g_yzz_x_y_xz, g_yzz_x_y_yy, g_yzz_x_y_yz, g_yzz_x_y_zz, g_zz_0_0_0_y_x_y_xx, g_zz_0_0_0_y_x_y_xy, g_zz_0_0_0_y_x_y_xz, g_zz_0_0_0_y_x_y_yy, g_zz_0_0_0_y_x_y_yz, g_zz_0_0_0_y_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_y_xx[i] = -2.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_yzz_x_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_y_xy[i] = -2.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_yzz_x_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_y_xz[i] = -2.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_yzz_x_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_y_yy[i] = -2.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_yzz_x_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_y_yz[i] = -2.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_yzz_x_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_y_zz[i] = -2.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_yzz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz, g_yzz_x_z_xx, g_yzz_x_z_xy, g_yzz_x_z_xz, g_yzz_x_z_yy, g_yzz_x_z_yz, g_yzz_x_z_zz, g_zz_0_0_0_y_x_z_xx, g_zz_0_0_0_y_x_z_xy, g_zz_0_0_0_y_x_z_xz, g_zz_0_0_0_y_x_z_yy, g_zz_0_0_0_y_x_z_yz, g_zz_0_0_0_y_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_z_xx[i] = -2.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_yzz_x_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_z_xy[i] = -2.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_yzz_x_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_z_xz[i] = -2.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_yzz_x_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_z_yy[i] = -2.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_yzz_x_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_z_yz[i] = -2.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_yzz_x_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_z_zz[i] = -2.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_yzz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz, g_yzz_y_x_xx, g_yzz_y_x_xy, g_yzz_y_x_xz, g_yzz_y_x_yy, g_yzz_y_x_yz, g_yzz_y_x_zz, g_zz_0_0_0_y_y_x_xx, g_zz_0_0_0_y_y_x_xy, g_zz_0_0_0_y_y_x_xz, g_zz_0_0_0_y_y_x_yy, g_zz_0_0_0_y_y_x_yz, g_zz_0_0_0_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_x_xx[i] = -2.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_yzz_y_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_x_xy[i] = -2.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_yzz_y_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_x_xz[i] = -2.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_yzz_y_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_x_yy[i] = -2.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_yzz_y_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_x_yz[i] = -2.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_yzz_y_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_x_zz[i] = -2.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_yzz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz, g_yzz_y_y_xx, g_yzz_y_y_xy, g_yzz_y_y_xz, g_yzz_y_y_yy, g_yzz_y_y_yz, g_yzz_y_y_zz, g_zz_0_0_0_y_y_y_xx, g_zz_0_0_0_y_y_y_xy, g_zz_0_0_0_y_y_y_xz, g_zz_0_0_0_y_y_y_yy, g_zz_0_0_0_y_y_y_yz, g_zz_0_0_0_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_y_xx[i] = -2.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_yzz_y_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_y_xy[i] = -2.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_yzz_y_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_y_xz[i] = -2.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_yzz_y_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_y_yy[i] = -2.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_yzz_y_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_y_yz[i] = -2.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_yzz_y_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_y_zz[i] = -2.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_yzz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz, g_yzz_y_z_xx, g_yzz_y_z_xy, g_yzz_y_z_xz, g_yzz_y_z_yy, g_yzz_y_z_yz, g_yzz_y_z_zz, g_zz_0_0_0_y_y_z_xx, g_zz_0_0_0_y_y_z_xy, g_zz_0_0_0_y_y_z_xz, g_zz_0_0_0_y_y_z_yy, g_zz_0_0_0_y_y_z_yz, g_zz_0_0_0_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_z_xx[i] = -2.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_yzz_y_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_z_xy[i] = -2.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_yzz_y_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_z_xz[i] = -2.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_yzz_y_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_z_yy[i] = -2.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_yzz_y_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_z_yz[i] = -2.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_yzz_y_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_z_zz[i] = -2.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_yzz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz, g_yzz_z_x_xx, g_yzz_z_x_xy, g_yzz_z_x_xz, g_yzz_z_x_yy, g_yzz_z_x_yz, g_yzz_z_x_zz, g_zz_0_0_0_y_z_x_xx, g_zz_0_0_0_y_z_x_xy, g_zz_0_0_0_y_z_x_xz, g_zz_0_0_0_y_z_x_yy, g_zz_0_0_0_y_z_x_yz, g_zz_0_0_0_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_x_xx[i] = -2.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_yzz_z_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_x_xy[i] = -2.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_yzz_z_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_x_xz[i] = -2.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_yzz_z_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_x_yy[i] = -2.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_yzz_z_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_x_yz[i] = -2.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_yzz_z_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_x_zz[i] = -2.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_yzz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz, g_yzz_z_y_xx, g_yzz_z_y_xy, g_yzz_z_y_xz, g_yzz_z_y_yy, g_yzz_z_y_yz, g_yzz_z_y_zz, g_zz_0_0_0_y_z_y_xx, g_zz_0_0_0_y_z_y_xy, g_zz_0_0_0_y_z_y_xz, g_zz_0_0_0_y_z_y_yy, g_zz_0_0_0_y_z_y_yz, g_zz_0_0_0_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_y_xx[i] = -2.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_yzz_z_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_y_xy[i] = -2.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_yzz_z_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_y_xz[i] = -2.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_yzz_z_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_y_yy[i] = -2.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_yzz_z_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_y_yz[i] = -2.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_yzz_z_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_y_zz[i] = -2.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_yzz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz, g_yzz_z_z_xx, g_yzz_z_z_xy, g_yzz_z_z_xz, g_yzz_z_z_yy, g_yzz_z_z_yz, g_yzz_z_z_zz, g_zz_0_0_0_y_z_z_xx, g_zz_0_0_0_y_z_z_xy, g_zz_0_0_0_y_z_z_xz, g_zz_0_0_0_y_z_z_yy, g_zz_0_0_0_y_z_z_yz, g_zz_0_0_0_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_z_xx[i] = -2.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_yzz_z_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_z_xy[i] = -2.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_yzz_z_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_z_xz[i] = -2.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_yzz_z_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_z_yy[i] = -2.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_yzz_z_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_z_yz[i] = -2.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_yzz_z_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_z_zz[i] = -2.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_yzz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz, g_zz_0_0_0_z_x_x_xx, g_zz_0_0_0_z_x_x_xy, g_zz_0_0_0_z_x_x_xz, g_zz_0_0_0_z_x_x_yy, g_zz_0_0_0_z_x_x_yz, g_zz_0_0_0_z_x_x_zz, g_zzz_x_x_xx, g_zzz_x_x_xy, g_zzz_x_x_xz, g_zzz_x_x_yy, g_zzz_x_x_yz, g_zzz_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_x_xx[i] = -6.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_zzz_x_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_x_xy[i] = -6.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_zzz_x_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_x_xz[i] = -6.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_zzz_x_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_x_yy[i] = -6.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_zzz_x_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_x_yz[i] = -6.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_zzz_x_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_x_zz[i] = -6.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_zzz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz, g_zz_0_0_0_z_x_y_xx, g_zz_0_0_0_z_x_y_xy, g_zz_0_0_0_z_x_y_xz, g_zz_0_0_0_z_x_y_yy, g_zz_0_0_0_z_x_y_yz, g_zz_0_0_0_z_x_y_zz, g_zzz_x_y_xx, g_zzz_x_y_xy, g_zzz_x_y_xz, g_zzz_x_y_yy, g_zzz_x_y_yz, g_zzz_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_y_xx[i] = -6.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_zzz_x_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_y_xy[i] = -6.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_zzz_x_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_y_xz[i] = -6.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_zzz_x_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_y_yy[i] = -6.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_zzz_x_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_y_yz[i] = -6.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_zzz_x_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_y_zz[i] = -6.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_zzz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz, g_zz_0_0_0_z_x_z_xx, g_zz_0_0_0_z_x_z_xy, g_zz_0_0_0_z_x_z_xz, g_zz_0_0_0_z_x_z_yy, g_zz_0_0_0_z_x_z_yz, g_zz_0_0_0_z_x_z_zz, g_zzz_x_z_xx, g_zzz_x_z_xy, g_zzz_x_z_xz, g_zzz_x_z_yy, g_zzz_x_z_yz, g_zzz_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_z_xx[i] = -6.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_zzz_x_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_z_xy[i] = -6.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_zzz_x_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_z_xz[i] = -6.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_zzz_x_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_z_yy[i] = -6.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_zzz_x_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_z_yz[i] = -6.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_zzz_x_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_z_zz[i] = -6.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_zzz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz, g_zz_0_0_0_z_y_x_xx, g_zz_0_0_0_z_y_x_xy, g_zz_0_0_0_z_y_x_xz, g_zz_0_0_0_z_y_x_yy, g_zz_0_0_0_z_y_x_yz, g_zz_0_0_0_z_y_x_zz, g_zzz_y_x_xx, g_zzz_y_x_xy, g_zzz_y_x_xz, g_zzz_y_x_yy, g_zzz_y_x_yz, g_zzz_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_x_xx[i] = -6.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_zzz_y_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_x_xy[i] = -6.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_zzz_y_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_x_xz[i] = -6.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_zzz_y_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_x_yy[i] = -6.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_zzz_y_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_x_yz[i] = -6.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_zzz_y_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_x_zz[i] = -6.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_zzz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz, g_zz_0_0_0_z_y_y_xx, g_zz_0_0_0_z_y_y_xy, g_zz_0_0_0_z_y_y_xz, g_zz_0_0_0_z_y_y_yy, g_zz_0_0_0_z_y_y_yz, g_zz_0_0_0_z_y_y_zz, g_zzz_y_y_xx, g_zzz_y_y_xy, g_zzz_y_y_xz, g_zzz_y_y_yy, g_zzz_y_y_yz, g_zzz_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_y_xx[i] = -6.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_zzz_y_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_y_xy[i] = -6.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_zzz_y_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_y_xz[i] = -6.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_zzz_y_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_y_yy[i] = -6.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_zzz_y_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_y_yz[i] = -6.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_zzz_y_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_y_zz[i] = -6.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_zzz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz, g_zz_0_0_0_z_y_z_xx, g_zz_0_0_0_z_y_z_xy, g_zz_0_0_0_z_y_z_xz, g_zz_0_0_0_z_y_z_yy, g_zz_0_0_0_z_y_z_yz, g_zz_0_0_0_z_y_z_zz, g_zzz_y_z_xx, g_zzz_y_z_xy, g_zzz_y_z_xz, g_zzz_y_z_yy, g_zzz_y_z_yz, g_zzz_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_z_xx[i] = -6.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_zzz_y_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_z_xy[i] = -6.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_zzz_y_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_z_xz[i] = -6.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_zzz_y_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_z_yy[i] = -6.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_zzz_y_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_z_yz[i] = -6.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_zzz_y_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_z_zz[i] = -6.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_zzz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz, g_zz_0_0_0_z_z_x_xx, g_zz_0_0_0_z_z_x_xy, g_zz_0_0_0_z_z_x_xz, g_zz_0_0_0_z_z_x_yy, g_zz_0_0_0_z_z_x_yz, g_zz_0_0_0_z_z_x_zz, g_zzz_z_x_xx, g_zzz_z_x_xy, g_zzz_z_x_xz, g_zzz_z_x_yy, g_zzz_z_x_yz, g_zzz_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_x_xx[i] = -6.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_zzz_z_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_x_xy[i] = -6.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_zzz_z_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_x_xz[i] = -6.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_zzz_z_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_x_yy[i] = -6.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_zzz_z_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_x_yz[i] = -6.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_zzz_z_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_x_zz[i] = -6.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_zzz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz, g_zz_0_0_0_z_z_y_xx, g_zz_0_0_0_z_z_y_xy, g_zz_0_0_0_z_z_y_xz, g_zz_0_0_0_z_z_y_yy, g_zz_0_0_0_z_z_y_yz, g_zz_0_0_0_z_z_y_zz, g_zzz_z_y_xx, g_zzz_z_y_xy, g_zzz_z_y_xz, g_zzz_z_y_yy, g_zzz_z_y_yz, g_zzz_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_y_xx[i] = -6.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_zzz_z_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_y_xy[i] = -6.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_zzz_z_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_y_xz[i] = -6.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_zzz_z_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_y_yy[i] = -6.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_zzz_z_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_y_yz[i] = -6.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_zzz_z_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_y_zz[i] = -6.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_zzz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz, g_zz_0_0_0_z_z_z_xx, g_zz_0_0_0_z_z_z_xy, g_zz_0_0_0_z_z_z_xz, g_zz_0_0_0_z_z_z_yy, g_zz_0_0_0_z_z_z_yz, g_zz_0_0_0_z_z_z_zz, g_zzz_z_z_xx, g_zzz_z_z_xy, g_zzz_z_z_xz, g_zzz_z_z_yy, g_zzz_z_z_yz, g_zzz_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_z_xx[i] = -6.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_zzz_z_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_z_xy[i] = -6.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_zzz_z_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_z_xz[i] = -6.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_zzz_z_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_z_yy[i] = -6.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_zzz_z_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_z_yz[i] = -6.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_zzz_z_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_z_zz[i] = -6.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_zzz_z_z_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

