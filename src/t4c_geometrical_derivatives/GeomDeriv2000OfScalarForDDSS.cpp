#include "GeomDeriv2000OfScalarForDDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_ddss_0(CSimdArray<double>& buffer_2000_ddss,
                     const CSimdArray<double>& buffer_sdss,
                     const CSimdArray<double>& buffer_ddss,
                     const CSimdArray<double>& buffer_gdss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_ddss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sdss

    auto g_0_xx_0_0 = buffer_sdss[0];

    auto g_0_xy_0_0 = buffer_sdss[1];

    auto g_0_xz_0_0 = buffer_sdss[2];

    auto g_0_yy_0_0 = buffer_sdss[3];

    auto g_0_yz_0_0 = buffer_sdss[4];

    auto g_0_zz_0_0 = buffer_sdss[5];

    /// Set up components of auxilary buffer : buffer_ddss

    auto g_xx_xx_0_0 = buffer_ddss[0];

    auto g_xx_xy_0_0 = buffer_ddss[1];

    auto g_xx_xz_0_0 = buffer_ddss[2];

    auto g_xx_yy_0_0 = buffer_ddss[3];

    auto g_xx_yz_0_0 = buffer_ddss[4];

    auto g_xx_zz_0_0 = buffer_ddss[5];

    auto g_xy_xx_0_0 = buffer_ddss[6];

    auto g_xy_xy_0_0 = buffer_ddss[7];

    auto g_xy_xz_0_0 = buffer_ddss[8];

    auto g_xy_yy_0_0 = buffer_ddss[9];

    auto g_xy_yz_0_0 = buffer_ddss[10];

    auto g_xy_zz_0_0 = buffer_ddss[11];

    auto g_xz_xx_0_0 = buffer_ddss[12];

    auto g_xz_xy_0_0 = buffer_ddss[13];

    auto g_xz_xz_0_0 = buffer_ddss[14];

    auto g_xz_yy_0_0 = buffer_ddss[15];

    auto g_xz_yz_0_0 = buffer_ddss[16];

    auto g_xz_zz_0_0 = buffer_ddss[17];

    auto g_yy_xx_0_0 = buffer_ddss[18];

    auto g_yy_xy_0_0 = buffer_ddss[19];

    auto g_yy_xz_0_0 = buffer_ddss[20];

    auto g_yy_yy_0_0 = buffer_ddss[21];

    auto g_yy_yz_0_0 = buffer_ddss[22];

    auto g_yy_zz_0_0 = buffer_ddss[23];

    auto g_yz_xx_0_0 = buffer_ddss[24];

    auto g_yz_xy_0_0 = buffer_ddss[25];

    auto g_yz_xz_0_0 = buffer_ddss[26];

    auto g_yz_yy_0_0 = buffer_ddss[27];

    auto g_yz_yz_0_0 = buffer_ddss[28];

    auto g_yz_zz_0_0 = buffer_ddss[29];

    auto g_zz_xx_0_0 = buffer_ddss[30];

    auto g_zz_xy_0_0 = buffer_ddss[31];

    auto g_zz_xz_0_0 = buffer_ddss[32];

    auto g_zz_yy_0_0 = buffer_ddss[33];

    auto g_zz_yz_0_0 = buffer_ddss[34];

    auto g_zz_zz_0_0 = buffer_ddss[35];

    /// Set up components of auxilary buffer : buffer_gdss

    auto g_xxxx_xx_0_0 = buffer_gdss[0];

    auto g_xxxx_xy_0_0 = buffer_gdss[1];

    auto g_xxxx_xz_0_0 = buffer_gdss[2];

    auto g_xxxx_yy_0_0 = buffer_gdss[3];

    auto g_xxxx_yz_0_0 = buffer_gdss[4];

    auto g_xxxx_zz_0_0 = buffer_gdss[5];

    auto g_xxxy_xx_0_0 = buffer_gdss[6];

    auto g_xxxy_xy_0_0 = buffer_gdss[7];

    auto g_xxxy_xz_0_0 = buffer_gdss[8];

    auto g_xxxy_yy_0_0 = buffer_gdss[9];

    auto g_xxxy_yz_0_0 = buffer_gdss[10];

    auto g_xxxy_zz_0_0 = buffer_gdss[11];

    auto g_xxxz_xx_0_0 = buffer_gdss[12];

    auto g_xxxz_xy_0_0 = buffer_gdss[13];

    auto g_xxxz_xz_0_0 = buffer_gdss[14];

    auto g_xxxz_yy_0_0 = buffer_gdss[15];

    auto g_xxxz_yz_0_0 = buffer_gdss[16];

    auto g_xxxz_zz_0_0 = buffer_gdss[17];

    auto g_xxyy_xx_0_0 = buffer_gdss[18];

    auto g_xxyy_xy_0_0 = buffer_gdss[19];

    auto g_xxyy_xz_0_0 = buffer_gdss[20];

    auto g_xxyy_yy_0_0 = buffer_gdss[21];

    auto g_xxyy_yz_0_0 = buffer_gdss[22];

    auto g_xxyy_zz_0_0 = buffer_gdss[23];

    auto g_xxyz_xx_0_0 = buffer_gdss[24];

    auto g_xxyz_xy_0_0 = buffer_gdss[25];

    auto g_xxyz_xz_0_0 = buffer_gdss[26];

    auto g_xxyz_yy_0_0 = buffer_gdss[27];

    auto g_xxyz_yz_0_0 = buffer_gdss[28];

    auto g_xxyz_zz_0_0 = buffer_gdss[29];

    auto g_xxzz_xx_0_0 = buffer_gdss[30];

    auto g_xxzz_xy_0_0 = buffer_gdss[31];

    auto g_xxzz_xz_0_0 = buffer_gdss[32];

    auto g_xxzz_yy_0_0 = buffer_gdss[33];

    auto g_xxzz_yz_0_0 = buffer_gdss[34];

    auto g_xxzz_zz_0_0 = buffer_gdss[35];

    auto g_xyyy_xx_0_0 = buffer_gdss[36];

    auto g_xyyy_xy_0_0 = buffer_gdss[37];

    auto g_xyyy_xz_0_0 = buffer_gdss[38];

    auto g_xyyy_yy_0_0 = buffer_gdss[39];

    auto g_xyyy_yz_0_0 = buffer_gdss[40];

    auto g_xyyy_zz_0_0 = buffer_gdss[41];

    auto g_xyyz_xx_0_0 = buffer_gdss[42];

    auto g_xyyz_xy_0_0 = buffer_gdss[43];

    auto g_xyyz_xz_0_0 = buffer_gdss[44];

    auto g_xyyz_yy_0_0 = buffer_gdss[45];

    auto g_xyyz_yz_0_0 = buffer_gdss[46];

    auto g_xyyz_zz_0_0 = buffer_gdss[47];

    auto g_xyzz_xx_0_0 = buffer_gdss[48];

    auto g_xyzz_xy_0_0 = buffer_gdss[49];

    auto g_xyzz_xz_0_0 = buffer_gdss[50];

    auto g_xyzz_yy_0_0 = buffer_gdss[51];

    auto g_xyzz_yz_0_0 = buffer_gdss[52];

    auto g_xyzz_zz_0_0 = buffer_gdss[53];

    auto g_xzzz_xx_0_0 = buffer_gdss[54];

    auto g_xzzz_xy_0_0 = buffer_gdss[55];

    auto g_xzzz_xz_0_0 = buffer_gdss[56];

    auto g_xzzz_yy_0_0 = buffer_gdss[57];

    auto g_xzzz_yz_0_0 = buffer_gdss[58];

    auto g_xzzz_zz_0_0 = buffer_gdss[59];

    auto g_yyyy_xx_0_0 = buffer_gdss[60];

    auto g_yyyy_xy_0_0 = buffer_gdss[61];

    auto g_yyyy_xz_0_0 = buffer_gdss[62];

    auto g_yyyy_yy_0_0 = buffer_gdss[63];

    auto g_yyyy_yz_0_0 = buffer_gdss[64];

    auto g_yyyy_zz_0_0 = buffer_gdss[65];

    auto g_yyyz_xx_0_0 = buffer_gdss[66];

    auto g_yyyz_xy_0_0 = buffer_gdss[67];

    auto g_yyyz_xz_0_0 = buffer_gdss[68];

    auto g_yyyz_yy_0_0 = buffer_gdss[69];

    auto g_yyyz_yz_0_0 = buffer_gdss[70];

    auto g_yyyz_zz_0_0 = buffer_gdss[71];

    auto g_yyzz_xx_0_0 = buffer_gdss[72];

    auto g_yyzz_xy_0_0 = buffer_gdss[73];

    auto g_yyzz_xz_0_0 = buffer_gdss[74];

    auto g_yyzz_yy_0_0 = buffer_gdss[75];

    auto g_yyzz_yz_0_0 = buffer_gdss[76];

    auto g_yyzz_zz_0_0 = buffer_gdss[77];

    auto g_yzzz_xx_0_0 = buffer_gdss[78];

    auto g_yzzz_xy_0_0 = buffer_gdss[79];

    auto g_yzzz_xz_0_0 = buffer_gdss[80];

    auto g_yzzz_yy_0_0 = buffer_gdss[81];

    auto g_yzzz_yz_0_0 = buffer_gdss[82];

    auto g_yzzz_zz_0_0 = buffer_gdss[83];

    auto g_zzzz_xx_0_0 = buffer_gdss[84];

    auto g_zzzz_xy_0_0 = buffer_gdss[85];

    auto g_zzzz_xz_0_0 = buffer_gdss[86];

    auto g_zzzz_yy_0_0 = buffer_gdss[87];

    auto g_zzzz_yz_0_0 = buffer_gdss[88];

    auto g_zzzz_zz_0_0 = buffer_gdss[89];

    /// Set up components of integrals buffer : buffer_2000_ddss

    auto g_xx_0_0_0_xx_xx_0_0 = buffer_2000_ddss[0];

    auto g_xx_0_0_0_xx_xy_0_0 = buffer_2000_ddss[1];

    auto g_xx_0_0_0_xx_xz_0_0 = buffer_2000_ddss[2];

    auto g_xx_0_0_0_xx_yy_0_0 = buffer_2000_ddss[3];

    auto g_xx_0_0_0_xx_yz_0_0 = buffer_2000_ddss[4];

    auto g_xx_0_0_0_xx_zz_0_0 = buffer_2000_ddss[5];

    auto g_xx_0_0_0_xy_xx_0_0 = buffer_2000_ddss[6];

    auto g_xx_0_0_0_xy_xy_0_0 = buffer_2000_ddss[7];

    auto g_xx_0_0_0_xy_xz_0_0 = buffer_2000_ddss[8];

    auto g_xx_0_0_0_xy_yy_0_0 = buffer_2000_ddss[9];

    auto g_xx_0_0_0_xy_yz_0_0 = buffer_2000_ddss[10];

    auto g_xx_0_0_0_xy_zz_0_0 = buffer_2000_ddss[11];

    auto g_xx_0_0_0_xz_xx_0_0 = buffer_2000_ddss[12];

    auto g_xx_0_0_0_xz_xy_0_0 = buffer_2000_ddss[13];

    auto g_xx_0_0_0_xz_xz_0_0 = buffer_2000_ddss[14];

    auto g_xx_0_0_0_xz_yy_0_0 = buffer_2000_ddss[15];

    auto g_xx_0_0_0_xz_yz_0_0 = buffer_2000_ddss[16];

    auto g_xx_0_0_0_xz_zz_0_0 = buffer_2000_ddss[17];

    auto g_xx_0_0_0_yy_xx_0_0 = buffer_2000_ddss[18];

    auto g_xx_0_0_0_yy_xy_0_0 = buffer_2000_ddss[19];

    auto g_xx_0_0_0_yy_xz_0_0 = buffer_2000_ddss[20];

    auto g_xx_0_0_0_yy_yy_0_0 = buffer_2000_ddss[21];

    auto g_xx_0_0_0_yy_yz_0_0 = buffer_2000_ddss[22];

    auto g_xx_0_0_0_yy_zz_0_0 = buffer_2000_ddss[23];

    auto g_xx_0_0_0_yz_xx_0_0 = buffer_2000_ddss[24];

    auto g_xx_0_0_0_yz_xy_0_0 = buffer_2000_ddss[25];

    auto g_xx_0_0_0_yz_xz_0_0 = buffer_2000_ddss[26];

    auto g_xx_0_0_0_yz_yy_0_0 = buffer_2000_ddss[27];

    auto g_xx_0_0_0_yz_yz_0_0 = buffer_2000_ddss[28];

    auto g_xx_0_0_0_yz_zz_0_0 = buffer_2000_ddss[29];

    auto g_xx_0_0_0_zz_xx_0_0 = buffer_2000_ddss[30];

    auto g_xx_0_0_0_zz_xy_0_0 = buffer_2000_ddss[31];

    auto g_xx_0_0_0_zz_xz_0_0 = buffer_2000_ddss[32];

    auto g_xx_0_0_0_zz_yy_0_0 = buffer_2000_ddss[33];

    auto g_xx_0_0_0_zz_yz_0_0 = buffer_2000_ddss[34];

    auto g_xx_0_0_0_zz_zz_0_0 = buffer_2000_ddss[35];

    auto g_xy_0_0_0_xx_xx_0_0 = buffer_2000_ddss[36];

    auto g_xy_0_0_0_xx_xy_0_0 = buffer_2000_ddss[37];

    auto g_xy_0_0_0_xx_xz_0_0 = buffer_2000_ddss[38];

    auto g_xy_0_0_0_xx_yy_0_0 = buffer_2000_ddss[39];

    auto g_xy_0_0_0_xx_yz_0_0 = buffer_2000_ddss[40];

    auto g_xy_0_0_0_xx_zz_0_0 = buffer_2000_ddss[41];

    auto g_xy_0_0_0_xy_xx_0_0 = buffer_2000_ddss[42];

    auto g_xy_0_0_0_xy_xy_0_0 = buffer_2000_ddss[43];

    auto g_xy_0_0_0_xy_xz_0_0 = buffer_2000_ddss[44];

    auto g_xy_0_0_0_xy_yy_0_0 = buffer_2000_ddss[45];

    auto g_xy_0_0_0_xy_yz_0_0 = buffer_2000_ddss[46];

    auto g_xy_0_0_0_xy_zz_0_0 = buffer_2000_ddss[47];

    auto g_xy_0_0_0_xz_xx_0_0 = buffer_2000_ddss[48];

    auto g_xy_0_0_0_xz_xy_0_0 = buffer_2000_ddss[49];

    auto g_xy_0_0_0_xz_xz_0_0 = buffer_2000_ddss[50];

    auto g_xy_0_0_0_xz_yy_0_0 = buffer_2000_ddss[51];

    auto g_xy_0_0_0_xz_yz_0_0 = buffer_2000_ddss[52];

    auto g_xy_0_0_0_xz_zz_0_0 = buffer_2000_ddss[53];

    auto g_xy_0_0_0_yy_xx_0_0 = buffer_2000_ddss[54];

    auto g_xy_0_0_0_yy_xy_0_0 = buffer_2000_ddss[55];

    auto g_xy_0_0_0_yy_xz_0_0 = buffer_2000_ddss[56];

    auto g_xy_0_0_0_yy_yy_0_0 = buffer_2000_ddss[57];

    auto g_xy_0_0_0_yy_yz_0_0 = buffer_2000_ddss[58];

    auto g_xy_0_0_0_yy_zz_0_0 = buffer_2000_ddss[59];

    auto g_xy_0_0_0_yz_xx_0_0 = buffer_2000_ddss[60];

    auto g_xy_0_0_0_yz_xy_0_0 = buffer_2000_ddss[61];

    auto g_xy_0_0_0_yz_xz_0_0 = buffer_2000_ddss[62];

    auto g_xy_0_0_0_yz_yy_0_0 = buffer_2000_ddss[63];

    auto g_xy_0_0_0_yz_yz_0_0 = buffer_2000_ddss[64];

    auto g_xy_0_0_0_yz_zz_0_0 = buffer_2000_ddss[65];

    auto g_xy_0_0_0_zz_xx_0_0 = buffer_2000_ddss[66];

    auto g_xy_0_0_0_zz_xy_0_0 = buffer_2000_ddss[67];

    auto g_xy_0_0_0_zz_xz_0_0 = buffer_2000_ddss[68];

    auto g_xy_0_0_0_zz_yy_0_0 = buffer_2000_ddss[69];

    auto g_xy_0_0_0_zz_yz_0_0 = buffer_2000_ddss[70];

    auto g_xy_0_0_0_zz_zz_0_0 = buffer_2000_ddss[71];

    auto g_xz_0_0_0_xx_xx_0_0 = buffer_2000_ddss[72];

    auto g_xz_0_0_0_xx_xy_0_0 = buffer_2000_ddss[73];

    auto g_xz_0_0_0_xx_xz_0_0 = buffer_2000_ddss[74];

    auto g_xz_0_0_0_xx_yy_0_0 = buffer_2000_ddss[75];

    auto g_xz_0_0_0_xx_yz_0_0 = buffer_2000_ddss[76];

    auto g_xz_0_0_0_xx_zz_0_0 = buffer_2000_ddss[77];

    auto g_xz_0_0_0_xy_xx_0_0 = buffer_2000_ddss[78];

    auto g_xz_0_0_0_xy_xy_0_0 = buffer_2000_ddss[79];

    auto g_xz_0_0_0_xy_xz_0_0 = buffer_2000_ddss[80];

    auto g_xz_0_0_0_xy_yy_0_0 = buffer_2000_ddss[81];

    auto g_xz_0_0_0_xy_yz_0_0 = buffer_2000_ddss[82];

    auto g_xz_0_0_0_xy_zz_0_0 = buffer_2000_ddss[83];

    auto g_xz_0_0_0_xz_xx_0_0 = buffer_2000_ddss[84];

    auto g_xz_0_0_0_xz_xy_0_0 = buffer_2000_ddss[85];

    auto g_xz_0_0_0_xz_xz_0_0 = buffer_2000_ddss[86];

    auto g_xz_0_0_0_xz_yy_0_0 = buffer_2000_ddss[87];

    auto g_xz_0_0_0_xz_yz_0_0 = buffer_2000_ddss[88];

    auto g_xz_0_0_0_xz_zz_0_0 = buffer_2000_ddss[89];

    auto g_xz_0_0_0_yy_xx_0_0 = buffer_2000_ddss[90];

    auto g_xz_0_0_0_yy_xy_0_0 = buffer_2000_ddss[91];

    auto g_xz_0_0_0_yy_xz_0_0 = buffer_2000_ddss[92];

    auto g_xz_0_0_0_yy_yy_0_0 = buffer_2000_ddss[93];

    auto g_xz_0_0_0_yy_yz_0_0 = buffer_2000_ddss[94];

    auto g_xz_0_0_0_yy_zz_0_0 = buffer_2000_ddss[95];

    auto g_xz_0_0_0_yz_xx_0_0 = buffer_2000_ddss[96];

    auto g_xz_0_0_0_yz_xy_0_0 = buffer_2000_ddss[97];

    auto g_xz_0_0_0_yz_xz_0_0 = buffer_2000_ddss[98];

    auto g_xz_0_0_0_yz_yy_0_0 = buffer_2000_ddss[99];

    auto g_xz_0_0_0_yz_yz_0_0 = buffer_2000_ddss[100];

    auto g_xz_0_0_0_yz_zz_0_0 = buffer_2000_ddss[101];

    auto g_xz_0_0_0_zz_xx_0_0 = buffer_2000_ddss[102];

    auto g_xz_0_0_0_zz_xy_0_0 = buffer_2000_ddss[103];

    auto g_xz_0_0_0_zz_xz_0_0 = buffer_2000_ddss[104];

    auto g_xz_0_0_0_zz_yy_0_0 = buffer_2000_ddss[105];

    auto g_xz_0_0_0_zz_yz_0_0 = buffer_2000_ddss[106];

    auto g_xz_0_0_0_zz_zz_0_0 = buffer_2000_ddss[107];

    auto g_yy_0_0_0_xx_xx_0_0 = buffer_2000_ddss[108];

    auto g_yy_0_0_0_xx_xy_0_0 = buffer_2000_ddss[109];

    auto g_yy_0_0_0_xx_xz_0_0 = buffer_2000_ddss[110];

    auto g_yy_0_0_0_xx_yy_0_0 = buffer_2000_ddss[111];

    auto g_yy_0_0_0_xx_yz_0_0 = buffer_2000_ddss[112];

    auto g_yy_0_0_0_xx_zz_0_0 = buffer_2000_ddss[113];

    auto g_yy_0_0_0_xy_xx_0_0 = buffer_2000_ddss[114];

    auto g_yy_0_0_0_xy_xy_0_0 = buffer_2000_ddss[115];

    auto g_yy_0_0_0_xy_xz_0_0 = buffer_2000_ddss[116];

    auto g_yy_0_0_0_xy_yy_0_0 = buffer_2000_ddss[117];

    auto g_yy_0_0_0_xy_yz_0_0 = buffer_2000_ddss[118];

    auto g_yy_0_0_0_xy_zz_0_0 = buffer_2000_ddss[119];

    auto g_yy_0_0_0_xz_xx_0_0 = buffer_2000_ddss[120];

    auto g_yy_0_0_0_xz_xy_0_0 = buffer_2000_ddss[121];

    auto g_yy_0_0_0_xz_xz_0_0 = buffer_2000_ddss[122];

    auto g_yy_0_0_0_xz_yy_0_0 = buffer_2000_ddss[123];

    auto g_yy_0_0_0_xz_yz_0_0 = buffer_2000_ddss[124];

    auto g_yy_0_0_0_xz_zz_0_0 = buffer_2000_ddss[125];

    auto g_yy_0_0_0_yy_xx_0_0 = buffer_2000_ddss[126];

    auto g_yy_0_0_0_yy_xy_0_0 = buffer_2000_ddss[127];

    auto g_yy_0_0_0_yy_xz_0_0 = buffer_2000_ddss[128];

    auto g_yy_0_0_0_yy_yy_0_0 = buffer_2000_ddss[129];

    auto g_yy_0_0_0_yy_yz_0_0 = buffer_2000_ddss[130];

    auto g_yy_0_0_0_yy_zz_0_0 = buffer_2000_ddss[131];

    auto g_yy_0_0_0_yz_xx_0_0 = buffer_2000_ddss[132];

    auto g_yy_0_0_0_yz_xy_0_0 = buffer_2000_ddss[133];

    auto g_yy_0_0_0_yz_xz_0_0 = buffer_2000_ddss[134];

    auto g_yy_0_0_0_yz_yy_0_0 = buffer_2000_ddss[135];

    auto g_yy_0_0_0_yz_yz_0_0 = buffer_2000_ddss[136];

    auto g_yy_0_0_0_yz_zz_0_0 = buffer_2000_ddss[137];

    auto g_yy_0_0_0_zz_xx_0_0 = buffer_2000_ddss[138];

    auto g_yy_0_0_0_zz_xy_0_0 = buffer_2000_ddss[139];

    auto g_yy_0_0_0_zz_xz_0_0 = buffer_2000_ddss[140];

    auto g_yy_0_0_0_zz_yy_0_0 = buffer_2000_ddss[141];

    auto g_yy_0_0_0_zz_yz_0_0 = buffer_2000_ddss[142];

    auto g_yy_0_0_0_zz_zz_0_0 = buffer_2000_ddss[143];

    auto g_yz_0_0_0_xx_xx_0_0 = buffer_2000_ddss[144];

    auto g_yz_0_0_0_xx_xy_0_0 = buffer_2000_ddss[145];

    auto g_yz_0_0_0_xx_xz_0_0 = buffer_2000_ddss[146];

    auto g_yz_0_0_0_xx_yy_0_0 = buffer_2000_ddss[147];

    auto g_yz_0_0_0_xx_yz_0_0 = buffer_2000_ddss[148];

    auto g_yz_0_0_0_xx_zz_0_0 = buffer_2000_ddss[149];

    auto g_yz_0_0_0_xy_xx_0_0 = buffer_2000_ddss[150];

    auto g_yz_0_0_0_xy_xy_0_0 = buffer_2000_ddss[151];

    auto g_yz_0_0_0_xy_xz_0_0 = buffer_2000_ddss[152];

    auto g_yz_0_0_0_xy_yy_0_0 = buffer_2000_ddss[153];

    auto g_yz_0_0_0_xy_yz_0_0 = buffer_2000_ddss[154];

    auto g_yz_0_0_0_xy_zz_0_0 = buffer_2000_ddss[155];

    auto g_yz_0_0_0_xz_xx_0_0 = buffer_2000_ddss[156];

    auto g_yz_0_0_0_xz_xy_0_0 = buffer_2000_ddss[157];

    auto g_yz_0_0_0_xz_xz_0_0 = buffer_2000_ddss[158];

    auto g_yz_0_0_0_xz_yy_0_0 = buffer_2000_ddss[159];

    auto g_yz_0_0_0_xz_yz_0_0 = buffer_2000_ddss[160];

    auto g_yz_0_0_0_xz_zz_0_0 = buffer_2000_ddss[161];

    auto g_yz_0_0_0_yy_xx_0_0 = buffer_2000_ddss[162];

    auto g_yz_0_0_0_yy_xy_0_0 = buffer_2000_ddss[163];

    auto g_yz_0_0_0_yy_xz_0_0 = buffer_2000_ddss[164];

    auto g_yz_0_0_0_yy_yy_0_0 = buffer_2000_ddss[165];

    auto g_yz_0_0_0_yy_yz_0_0 = buffer_2000_ddss[166];

    auto g_yz_0_0_0_yy_zz_0_0 = buffer_2000_ddss[167];

    auto g_yz_0_0_0_yz_xx_0_0 = buffer_2000_ddss[168];

    auto g_yz_0_0_0_yz_xy_0_0 = buffer_2000_ddss[169];

    auto g_yz_0_0_0_yz_xz_0_0 = buffer_2000_ddss[170];

    auto g_yz_0_0_0_yz_yy_0_0 = buffer_2000_ddss[171];

    auto g_yz_0_0_0_yz_yz_0_0 = buffer_2000_ddss[172];

    auto g_yz_0_0_0_yz_zz_0_0 = buffer_2000_ddss[173];

    auto g_yz_0_0_0_zz_xx_0_0 = buffer_2000_ddss[174];

    auto g_yz_0_0_0_zz_xy_0_0 = buffer_2000_ddss[175];

    auto g_yz_0_0_0_zz_xz_0_0 = buffer_2000_ddss[176];

    auto g_yz_0_0_0_zz_yy_0_0 = buffer_2000_ddss[177];

    auto g_yz_0_0_0_zz_yz_0_0 = buffer_2000_ddss[178];

    auto g_yz_0_0_0_zz_zz_0_0 = buffer_2000_ddss[179];

    auto g_zz_0_0_0_xx_xx_0_0 = buffer_2000_ddss[180];

    auto g_zz_0_0_0_xx_xy_0_0 = buffer_2000_ddss[181];

    auto g_zz_0_0_0_xx_xz_0_0 = buffer_2000_ddss[182];

    auto g_zz_0_0_0_xx_yy_0_0 = buffer_2000_ddss[183];

    auto g_zz_0_0_0_xx_yz_0_0 = buffer_2000_ddss[184];

    auto g_zz_0_0_0_xx_zz_0_0 = buffer_2000_ddss[185];

    auto g_zz_0_0_0_xy_xx_0_0 = buffer_2000_ddss[186];

    auto g_zz_0_0_0_xy_xy_0_0 = buffer_2000_ddss[187];

    auto g_zz_0_0_0_xy_xz_0_0 = buffer_2000_ddss[188];

    auto g_zz_0_0_0_xy_yy_0_0 = buffer_2000_ddss[189];

    auto g_zz_0_0_0_xy_yz_0_0 = buffer_2000_ddss[190];

    auto g_zz_0_0_0_xy_zz_0_0 = buffer_2000_ddss[191];

    auto g_zz_0_0_0_xz_xx_0_0 = buffer_2000_ddss[192];

    auto g_zz_0_0_0_xz_xy_0_0 = buffer_2000_ddss[193];

    auto g_zz_0_0_0_xz_xz_0_0 = buffer_2000_ddss[194];

    auto g_zz_0_0_0_xz_yy_0_0 = buffer_2000_ddss[195];

    auto g_zz_0_0_0_xz_yz_0_0 = buffer_2000_ddss[196];

    auto g_zz_0_0_0_xz_zz_0_0 = buffer_2000_ddss[197];

    auto g_zz_0_0_0_yy_xx_0_0 = buffer_2000_ddss[198];

    auto g_zz_0_0_0_yy_xy_0_0 = buffer_2000_ddss[199];

    auto g_zz_0_0_0_yy_xz_0_0 = buffer_2000_ddss[200];

    auto g_zz_0_0_0_yy_yy_0_0 = buffer_2000_ddss[201];

    auto g_zz_0_0_0_yy_yz_0_0 = buffer_2000_ddss[202];

    auto g_zz_0_0_0_yy_zz_0_0 = buffer_2000_ddss[203];

    auto g_zz_0_0_0_yz_xx_0_0 = buffer_2000_ddss[204];

    auto g_zz_0_0_0_yz_xy_0_0 = buffer_2000_ddss[205];

    auto g_zz_0_0_0_yz_xz_0_0 = buffer_2000_ddss[206];

    auto g_zz_0_0_0_yz_yy_0_0 = buffer_2000_ddss[207];

    auto g_zz_0_0_0_yz_yz_0_0 = buffer_2000_ddss[208];

    auto g_zz_0_0_0_yz_zz_0_0 = buffer_2000_ddss[209];

    auto g_zz_0_0_0_zz_xx_0_0 = buffer_2000_ddss[210];

    auto g_zz_0_0_0_zz_xy_0_0 = buffer_2000_ddss[211];

    auto g_zz_0_0_0_zz_xz_0_0 = buffer_2000_ddss[212];

    auto g_zz_0_0_0_zz_yy_0_0 = buffer_2000_ddss[213];

    auto g_zz_0_0_0_zz_yz_0_0 = buffer_2000_ddss[214];

    auto g_zz_0_0_0_zz_zz_0_0 = buffer_2000_ddss[215];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_xx_0_0_0_xx_xx_0_0, g_xx_0_0_0_xx_xy_0_0, g_xx_0_0_0_xx_xz_0_0, g_xx_0_0_0_xx_yy_0_0, g_xx_0_0_0_xx_yz_0_0, g_xx_0_0_0_xx_zz_0_0, g_xx_xx_0_0, g_xx_xy_0_0, g_xx_xz_0_0, g_xx_yy_0_0, g_xx_yz_0_0, g_xx_zz_0_0, g_xxxx_xx_0_0, g_xxxx_xy_0_0, g_xxxx_xz_0_0, g_xxxx_yy_0_0, g_xxxx_yz_0_0, g_xxxx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xx_0_0[i] = 2.0 * g_0_xx_0_0[i] - 10.0 * g_xx_xx_0_0[i] * a_exp + 4.0 * g_xxxx_xx_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_0_0[i] = 2.0 * g_0_xy_0_0[i] - 10.0 * g_xx_xy_0_0[i] * a_exp + 4.0 * g_xxxx_xy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_0_0[i] = 2.0 * g_0_xz_0_0[i] - 10.0 * g_xx_xz_0_0[i] * a_exp + 4.0 * g_xxxx_xz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_0_0[i] = 2.0 * g_0_yy_0_0[i] - 10.0 * g_xx_yy_0_0[i] * a_exp + 4.0 * g_xxxx_yy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_0_0[i] = 2.0 * g_0_yz_0_0[i] - 10.0 * g_xx_yz_0_0[i] * a_exp + 4.0 * g_xxxx_yz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_0_0[i] = 2.0 * g_0_zz_0_0[i] - 10.0 * g_xx_zz_0_0[i] * a_exp + 4.0 * g_xxxx_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xx_0_0, g_xx_0_0_0_xy_xy_0_0, g_xx_0_0_0_xy_xz_0_0, g_xx_0_0_0_xy_yy_0_0, g_xx_0_0_0_xy_yz_0_0, g_xx_0_0_0_xy_zz_0_0, g_xxxy_xx_0_0, g_xxxy_xy_0_0, g_xxxy_xz_0_0, g_xxxy_yy_0_0, g_xxxy_yz_0_0, g_xxxy_zz_0_0, g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xx_0_0[i] = -6.0 * g_xy_xx_0_0[i] * a_exp + 4.0 * g_xxxy_xx_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_0_0[i] = -6.0 * g_xy_xy_0_0[i] * a_exp + 4.0 * g_xxxy_xy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_0_0[i] = -6.0 * g_xy_xz_0_0[i] * a_exp + 4.0 * g_xxxy_xz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_0_0[i] = -6.0 * g_xy_yy_0_0[i] * a_exp + 4.0 * g_xxxy_yy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_0_0[i] = -6.0 * g_xy_yz_0_0[i] * a_exp + 4.0 * g_xxxy_yz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_0_0[i] = -6.0 * g_xy_zz_0_0[i] * a_exp + 4.0 * g_xxxy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xx_0_0, g_xx_0_0_0_xz_xy_0_0, g_xx_0_0_0_xz_xz_0_0, g_xx_0_0_0_xz_yy_0_0, g_xx_0_0_0_xz_yz_0_0, g_xx_0_0_0_xz_zz_0_0, g_xxxz_xx_0_0, g_xxxz_xy_0_0, g_xxxz_xz_0_0, g_xxxz_yy_0_0, g_xxxz_yz_0_0, g_xxxz_zz_0_0, g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xx_0_0[i] = -6.0 * g_xz_xx_0_0[i] * a_exp + 4.0 * g_xxxz_xx_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_0_0[i] = -6.0 * g_xz_xy_0_0[i] * a_exp + 4.0 * g_xxxz_xy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_0_0[i] = -6.0 * g_xz_xz_0_0[i] * a_exp + 4.0 * g_xxxz_xz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_0_0[i] = -6.0 * g_xz_yy_0_0[i] * a_exp + 4.0 * g_xxxz_yy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_0_0[i] = -6.0 * g_xz_yz_0_0[i] * a_exp + 4.0 * g_xxxz_yz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_0_0[i] = -6.0 * g_xz_zz_0_0[i] * a_exp + 4.0 * g_xxxz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xx_0_0, g_xx_0_0_0_yy_xy_0_0, g_xx_0_0_0_yy_xz_0_0, g_xx_0_0_0_yy_yy_0_0, g_xx_0_0_0_yy_yz_0_0, g_xx_0_0_0_yy_zz_0_0, g_xxyy_xx_0_0, g_xxyy_xy_0_0, g_xxyy_xz_0_0, g_xxyy_yy_0_0, g_xxyy_yz_0_0, g_xxyy_zz_0_0, g_yy_xx_0_0, g_yy_xy_0_0, g_yy_xz_0_0, g_yy_yy_0_0, g_yy_yz_0_0, g_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xx_0_0[i] = -2.0 * g_yy_xx_0_0[i] * a_exp + 4.0 * g_xxyy_xx_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_0_0[i] = -2.0 * g_yy_xy_0_0[i] * a_exp + 4.0 * g_xxyy_xy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_0_0[i] = -2.0 * g_yy_xz_0_0[i] * a_exp + 4.0 * g_xxyy_xz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_0_0[i] = -2.0 * g_yy_yy_0_0[i] * a_exp + 4.0 * g_xxyy_yy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_0_0[i] = -2.0 * g_yy_yz_0_0[i] * a_exp + 4.0 * g_xxyy_yz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_0_0[i] = -2.0 * g_yy_zz_0_0[i] * a_exp + 4.0 * g_xxyy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xx_0_0, g_xx_0_0_0_yz_xy_0_0, g_xx_0_0_0_yz_xz_0_0, g_xx_0_0_0_yz_yy_0_0, g_xx_0_0_0_yz_yz_0_0, g_xx_0_0_0_yz_zz_0_0, g_xxyz_xx_0_0, g_xxyz_xy_0_0, g_xxyz_xz_0_0, g_xxyz_yy_0_0, g_xxyz_yz_0_0, g_xxyz_zz_0_0, g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xx_0_0[i] = -2.0 * g_yz_xx_0_0[i] * a_exp + 4.0 * g_xxyz_xx_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_0_0[i] = -2.0 * g_yz_xy_0_0[i] * a_exp + 4.0 * g_xxyz_xy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_0_0[i] = -2.0 * g_yz_xz_0_0[i] * a_exp + 4.0 * g_xxyz_xz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_0_0[i] = -2.0 * g_yz_yy_0_0[i] * a_exp + 4.0 * g_xxyz_yy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_0_0[i] = -2.0 * g_yz_yz_0_0[i] * a_exp + 4.0 * g_xxyz_yz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_0_0[i] = -2.0 * g_yz_zz_0_0[i] * a_exp + 4.0 * g_xxyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xx_0_0, g_xx_0_0_0_zz_xy_0_0, g_xx_0_0_0_zz_xz_0_0, g_xx_0_0_0_zz_yy_0_0, g_xx_0_0_0_zz_yz_0_0, g_xx_0_0_0_zz_zz_0_0, g_xxzz_xx_0_0, g_xxzz_xy_0_0, g_xxzz_xz_0_0, g_xxzz_yy_0_0, g_xxzz_yz_0_0, g_xxzz_zz_0_0, g_zz_xx_0_0, g_zz_xy_0_0, g_zz_xz_0_0, g_zz_yy_0_0, g_zz_yz_0_0, g_zz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xx_0_0[i] = -2.0 * g_zz_xx_0_0[i] * a_exp + 4.0 * g_xxzz_xx_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_0_0[i] = -2.0 * g_zz_xy_0_0[i] * a_exp + 4.0 * g_xxzz_xy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_0_0[i] = -2.0 * g_zz_xz_0_0[i] * a_exp + 4.0 * g_xxzz_xz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_0_0[i] = -2.0 * g_zz_yy_0_0[i] * a_exp + 4.0 * g_xxzz_yy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_0_0[i] = -2.0 * g_zz_yz_0_0[i] * a_exp + 4.0 * g_xxzz_yz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_0_0[i] = -2.0 * g_zz_zz_0_0[i] * a_exp + 4.0 * g_xxzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xxxy_xx_0_0, g_xxxy_xy_0_0, g_xxxy_xz_0_0, g_xxxy_yy_0_0, g_xxxy_yz_0_0, g_xxxy_zz_0_0, g_xy_0_0_0_xx_xx_0_0, g_xy_0_0_0_xx_xy_0_0, g_xy_0_0_0_xx_xz_0_0, g_xy_0_0_0_xx_yy_0_0, g_xy_0_0_0_xx_yz_0_0, g_xy_0_0_0_xx_zz_0_0, g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xx_0_0[i] = -4.0 * g_xy_xx_0_0[i] * a_exp + 4.0 * g_xxxy_xx_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_0_0[i] = -4.0 * g_xy_xy_0_0[i] * a_exp + 4.0 * g_xxxy_xy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_0_0[i] = -4.0 * g_xy_xz_0_0[i] * a_exp + 4.0 * g_xxxy_xz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_0_0[i] = -4.0 * g_xy_yy_0_0[i] * a_exp + 4.0 * g_xxxy_yy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_0_0[i] = -4.0 * g_xy_yz_0_0[i] * a_exp + 4.0 * g_xxxy_yz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_0_0[i] = -4.0 * g_xy_zz_0_0[i] * a_exp + 4.0 * g_xxxy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_xx_xx_0_0, g_xx_xy_0_0, g_xx_xz_0_0, g_xx_yy_0_0, g_xx_yz_0_0, g_xx_zz_0_0, g_xxyy_xx_0_0, g_xxyy_xy_0_0, g_xxyy_xz_0_0, g_xxyy_yy_0_0, g_xxyy_yz_0_0, g_xxyy_zz_0_0, g_xy_0_0_0_xy_xx_0_0, g_xy_0_0_0_xy_xy_0_0, g_xy_0_0_0_xy_xz_0_0, g_xy_0_0_0_xy_yy_0_0, g_xy_0_0_0_xy_yz_0_0, g_xy_0_0_0_xy_zz_0_0, g_yy_xx_0_0, g_yy_xy_0_0, g_yy_xz_0_0, g_yy_yy_0_0, g_yy_yz_0_0, g_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xx_0_0[i] = g_0_xx_0_0[i] - 2.0 * g_yy_xx_0_0[i] * a_exp - 2.0 * g_xx_xx_0_0[i] * a_exp + 4.0 * g_xxyy_xx_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_0_0[i] = g_0_xy_0_0[i] - 2.0 * g_yy_xy_0_0[i] * a_exp - 2.0 * g_xx_xy_0_0[i] * a_exp + 4.0 * g_xxyy_xy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_0_0[i] = g_0_xz_0_0[i] - 2.0 * g_yy_xz_0_0[i] * a_exp - 2.0 * g_xx_xz_0_0[i] * a_exp + 4.0 * g_xxyy_xz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_0_0[i] = g_0_yy_0_0[i] - 2.0 * g_yy_yy_0_0[i] * a_exp - 2.0 * g_xx_yy_0_0[i] * a_exp + 4.0 * g_xxyy_yy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_0_0[i] = g_0_yz_0_0[i] - 2.0 * g_yy_yz_0_0[i] * a_exp - 2.0 * g_xx_yz_0_0[i] * a_exp + 4.0 * g_xxyy_yz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_0_0[i] = g_0_zz_0_0[i] - 2.0 * g_yy_zz_0_0[i] * a_exp - 2.0 * g_xx_zz_0_0[i] * a_exp + 4.0 * g_xxyy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_xxyz_xx_0_0, g_xxyz_xy_0_0, g_xxyz_xz_0_0, g_xxyz_yy_0_0, g_xxyz_yz_0_0, g_xxyz_zz_0_0, g_xy_0_0_0_xz_xx_0_0, g_xy_0_0_0_xz_xy_0_0, g_xy_0_0_0_xz_xz_0_0, g_xy_0_0_0_xz_yy_0_0, g_xy_0_0_0_xz_yz_0_0, g_xy_0_0_0_xz_zz_0_0, g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xx_0_0[i] = -2.0 * g_yz_xx_0_0[i] * a_exp + 4.0 * g_xxyz_xx_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_0_0[i] = -2.0 * g_yz_xy_0_0[i] * a_exp + 4.0 * g_xxyz_xy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_0_0[i] = -2.0 * g_yz_xz_0_0[i] * a_exp + 4.0 * g_xxyz_xz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_0_0[i] = -2.0 * g_yz_yy_0_0[i] * a_exp + 4.0 * g_xxyz_yy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_0_0[i] = -2.0 * g_yz_yz_0_0[i] * a_exp + 4.0 * g_xxyz_yz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_0_0[i] = -2.0 * g_yz_zz_0_0[i] * a_exp + 4.0 * g_xxyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xx_0_0, g_xy_0_0_0_yy_xy_0_0, g_xy_0_0_0_yy_xz_0_0, g_xy_0_0_0_yy_yy_0_0, g_xy_0_0_0_yy_yz_0_0, g_xy_0_0_0_yy_zz_0_0, g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_xy_zz_0_0, g_xyyy_xx_0_0, g_xyyy_xy_0_0, g_xyyy_xz_0_0, g_xyyy_yy_0_0, g_xyyy_yz_0_0, g_xyyy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xx_0_0[i] = -4.0 * g_xy_xx_0_0[i] * a_exp + 4.0 * g_xyyy_xx_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_0_0[i] = -4.0 * g_xy_xy_0_0[i] * a_exp + 4.0 * g_xyyy_xy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_0_0[i] = -4.0 * g_xy_xz_0_0[i] * a_exp + 4.0 * g_xyyy_xz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_0_0[i] = -4.0 * g_xy_yy_0_0[i] * a_exp + 4.0 * g_xyyy_yy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_0_0[i] = -4.0 * g_xy_yz_0_0[i] * a_exp + 4.0 * g_xyyy_yz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_0_0[i] = -4.0 * g_xy_zz_0_0[i] * a_exp + 4.0 * g_xyyy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xx_0_0, g_xy_0_0_0_yz_xy_0_0, g_xy_0_0_0_yz_xz_0_0, g_xy_0_0_0_yz_yy_0_0, g_xy_0_0_0_yz_yz_0_0, g_xy_0_0_0_yz_zz_0_0, g_xyyz_xx_0_0, g_xyyz_xy_0_0, g_xyyz_xz_0_0, g_xyyz_yy_0_0, g_xyyz_yz_0_0, g_xyyz_zz_0_0, g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xx_0_0[i] = -2.0 * g_xz_xx_0_0[i] * a_exp + 4.0 * g_xyyz_xx_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_0_0[i] = -2.0 * g_xz_xy_0_0[i] * a_exp + 4.0 * g_xyyz_xy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_0_0[i] = -2.0 * g_xz_xz_0_0[i] * a_exp + 4.0 * g_xyyz_xz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_0_0[i] = -2.0 * g_xz_yy_0_0[i] * a_exp + 4.0 * g_xyyz_yy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_0_0[i] = -2.0 * g_xz_yz_0_0[i] * a_exp + 4.0 * g_xyyz_yz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_0_0[i] = -2.0 * g_xz_zz_0_0[i] * a_exp + 4.0 * g_xyyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xx_0_0, g_xy_0_0_0_zz_xy_0_0, g_xy_0_0_0_zz_xz_0_0, g_xy_0_0_0_zz_yy_0_0, g_xy_0_0_0_zz_yz_0_0, g_xy_0_0_0_zz_zz_0_0, g_xyzz_xx_0_0, g_xyzz_xy_0_0, g_xyzz_xz_0_0, g_xyzz_yy_0_0, g_xyzz_yz_0_0, g_xyzz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xx_0_0[i] = 4.0 * g_xyzz_xx_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_0_0[i] = 4.0 * g_xyzz_xy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_0_0[i] = 4.0 * g_xyzz_xz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_0_0[i] = 4.0 * g_xyzz_yy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_0_0[i] = 4.0 * g_xyzz_yz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_0_0[i] = 4.0 * g_xyzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xxxz_xx_0_0, g_xxxz_xy_0_0, g_xxxz_xz_0_0, g_xxxz_yy_0_0, g_xxxz_yz_0_0, g_xxxz_zz_0_0, g_xz_0_0_0_xx_xx_0_0, g_xz_0_0_0_xx_xy_0_0, g_xz_0_0_0_xx_xz_0_0, g_xz_0_0_0_xx_yy_0_0, g_xz_0_0_0_xx_yz_0_0, g_xz_0_0_0_xx_zz_0_0, g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xx_0_0[i] = -4.0 * g_xz_xx_0_0[i] * a_exp + 4.0 * g_xxxz_xx_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_0_0[i] = -4.0 * g_xz_xy_0_0[i] * a_exp + 4.0 * g_xxxz_xy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_0_0[i] = -4.0 * g_xz_xz_0_0[i] * a_exp + 4.0 * g_xxxz_xz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_0_0[i] = -4.0 * g_xz_yy_0_0[i] * a_exp + 4.0 * g_xxxz_yy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_0_0[i] = -4.0 * g_xz_yz_0_0[i] * a_exp + 4.0 * g_xxxz_yz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_0_0[i] = -4.0 * g_xz_zz_0_0[i] * a_exp + 4.0 * g_xxxz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_xxyz_xx_0_0, g_xxyz_xy_0_0, g_xxyz_xz_0_0, g_xxyz_yy_0_0, g_xxyz_yz_0_0, g_xxyz_zz_0_0, g_xz_0_0_0_xy_xx_0_0, g_xz_0_0_0_xy_xy_0_0, g_xz_0_0_0_xy_xz_0_0, g_xz_0_0_0_xy_yy_0_0, g_xz_0_0_0_xy_yz_0_0, g_xz_0_0_0_xy_zz_0_0, g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xx_0_0[i] = -2.0 * g_yz_xx_0_0[i] * a_exp + 4.0 * g_xxyz_xx_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_0_0[i] = -2.0 * g_yz_xy_0_0[i] * a_exp + 4.0 * g_xxyz_xy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_0_0[i] = -2.0 * g_yz_xz_0_0[i] * a_exp + 4.0 * g_xxyz_xz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_0_0[i] = -2.0 * g_yz_yy_0_0[i] * a_exp + 4.0 * g_xxyz_yy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_0_0[i] = -2.0 * g_yz_yz_0_0[i] * a_exp + 4.0 * g_xxyz_yz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_0_0[i] = -2.0 * g_yz_zz_0_0[i] * a_exp + 4.0 * g_xxyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_xx_xx_0_0, g_xx_xy_0_0, g_xx_xz_0_0, g_xx_yy_0_0, g_xx_yz_0_0, g_xx_zz_0_0, g_xxzz_xx_0_0, g_xxzz_xy_0_0, g_xxzz_xz_0_0, g_xxzz_yy_0_0, g_xxzz_yz_0_0, g_xxzz_zz_0_0, g_xz_0_0_0_xz_xx_0_0, g_xz_0_0_0_xz_xy_0_0, g_xz_0_0_0_xz_xz_0_0, g_xz_0_0_0_xz_yy_0_0, g_xz_0_0_0_xz_yz_0_0, g_xz_0_0_0_xz_zz_0_0, g_zz_xx_0_0, g_zz_xy_0_0, g_zz_xz_0_0, g_zz_yy_0_0, g_zz_yz_0_0, g_zz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xx_0_0[i] = g_0_xx_0_0[i] - 2.0 * g_zz_xx_0_0[i] * a_exp - 2.0 * g_xx_xx_0_0[i] * a_exp + 4.0 * g_xxzz_xx_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_0_0[i] = g_0_xy_0_0[i] - 2.0 * g_zz_xy_0_0[i] * a_exp - 2.0 * g_xx_xy_0_0[i] * a_exp + 4.0 * g_xxzz_xy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_0_0[i] = g_0_xz_0_0[i] - 2.0 * g_zz_xz_0_0[i] * a_exp - 2.0 * g_xx_xz_0_0[i] * a_exp + 4.0 * g_xxzz_xz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_0_0[i] = g_0_yy_0_0[i] - 2.0 * g_zz_yy_0_0[i] * a_exp - 2.0 * g_xx_yy_0_0[i] * a_exp + 4.0 * g_xxzz_yy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_0_0[i] = g_0_yz_0_0[i] - 2.0 * g_zz_yz_0_0[i] * a_exp - 2.0 * g_xx_yz_0_0[i] * a_exp + 4.0 * g_xxzz_yz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_0_0[i] = g_0_zz_0_0[i] - 2.0 * g_zz_zz_0_0[i] * a_exp - 2.0 * g_xx_zz_0_0[i] * a_exp + 4.0 * g_xxzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_xyyz_xx_0_0, g_xyyz_xy_0_0, g_xyyz_xz_0_0, g_xyyz_yy_0_0, g_xyyz_yz_0_0, g_xyyz_zz_0_0, g_xz_0_0_0_yy_xx_0_0, g_xz_0_0_0_yy_xy_0_0, g_xz_0_0_0_yy_xz_0_0, g_xz_0_0_0_yy_yy_0_0, g_xz_0_0_0_yy_yz_0_0, g_xz_0_0_0_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xx_0_0[i] = 4.0 * g_xyyz_xx_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_0_0[i] = 4.0 * g_xyyz_xy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_0_0[i] = 4.0 * g_xyyz_xz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_0_0[i] = 4.0 * g_xyyz_yy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_0_0[i] = 4.0 * g_xyyz_yz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_0_0[i] = 4.0 * g_xyyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_xy_zz_0_0, g_xyzz_xx_0_0, g_xyzz_xy_0_0, g_xyzz_xz_0_0, g_xyzz_yy_0_0, g_xyzz_yz_0_0, g_xyzz_zz_0_0, g_xz_0_0_0_yz_xx_0_0, g_xz_0_0_0_yz_xy_0_0, g_xz_0_0_0_yz_xz_0_0, g_xz_0_0_0_yz_yy_0_0, g_xz_0_0_0_yz_yz_0_0, g_xz_0_0_0_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xx_0_0[i] = -2.0 * g_xy_xx_0_0[i] * a_exp + 4.0 * g_xyzz_xx_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_0_0[i] = -2.0 * g_xy_xy_0_0[i] * a_exp + 4.0 * g_xyzz_xy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_0_0[i] = -2.0 * g_xy_xz_0_0[i] * a_exp + 4.0 * g_xyzz_xz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_0_0[i] = -2.0 * g_xy_yy_0_0[i] * a_exp + 4.0 * g_xyzz_yy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_0_0[i] = -2.0 * g_xy_yz_0_0[i] * a_exp + 4.0 * g_xyzz_yz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_0_0[i] = -2.0 * g_xy_zz_0_0[i] * a_exp + 4.0 * g_xyzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xx_0_0, g_xz_0_0_0_zz_xy_0_0, g_xz_0_0_0_zz_xz_0_0, g_xz_0_0_0_zz_yy_0_0, g_xz_0_0_0_zz_yz_0_0, g_xz_0_0_0_zz_zz_0_0, g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_xz_zz_0_0, g_xzzz_xx_0_0, g_xzzz_xy_0_0, g_xzzz_xz_0_0, g_xzzz_yy_0_0, g_xzzz_yz_0_0, g_xzzz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xx_0_0[i] = -4.0 * g_xz_xx_0_0[i] * a_exp + 4.0 * g_xzzz_xx_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_0_0[i] = -4.0 * g_xz_xy_0_0[i] * a_exp + 4.0 * g_xzzz_xy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_0_0[i] = -4.0 * g_xz_xz_0_0[i] * a_exp + 4.0 * g_xzzz_xz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_0_0[i] = -4.0 * g_xz_yy_0_0[i] * a_exp + 4.0 * g_xzzz_yy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_0_0[i] = -4.0 * g_xz_yz_0_0[i] * a_exp + 4.0 * g_xzzz_yz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_0_0[i] = -4.0 * g_xz_zz_0_0[i] * a_exp + 4.0 * g_xzzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xx_xx_0_0, g_xx_xy_0_0, g_xx_xz_0_0, g_xx_yy_0_0, g_xx_yz_0_0, g_xx_zz_0_0, g_xxyy_xx_0_0, g_xxyy_xy_0_0, g_xxyy_xz_0_0, g_xxyy_yy_0_0, g_xxyy_yz_0_0, g_xxyy_zz_0_0, g_yy_0_0_0_xx_xx_0_0, g_yy_0_0_0_xx_xy_0_0, g_yy_0_0_0_xx_xz_0_0, g_yy_0_0_0_xx_yy_0_0, g_yy_0_0_0_xx_yz_0_0, g_yy_0_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xx_0_0[i] = -2.0 * g_xx_xx_0_0[i] * a_exp + 4.0 * g_xxyy_xx_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_0_0[i] = -2.0 * g_xx_xy_0_0[i] * a_exp + 4.0 * g_xxyy_xy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_0_0[i] = -2.0 * g_xx_xz_0_0[i] * a_exp + 4.0 * g_xxyy_xz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_0_0[i] = -2.0 * g_xx_yy_0_0[i] * a_exp + 4.0 * g_xxyy_yy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_0_0[i] = -2.0 * g_xx_yz_0_0[i] * a_exp + 4.0 * g_xxyy_yz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_0_0[i] = -2.0 * g_xx_zz_0_0[i] * a_exp + 4.0 * g_xxyy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_xy_zz_0_0, g_xyyy_xx_0_0, g_xyyy_xy_0_0, g_xyyy_xz_0_0, g_xyyy_yy_0_0, g_xyyy_yz_0_0, g_xyyy_zz_0_0, g_yy_0_0_0_xy_xx_0_0, g_yy_0_0_0_xy_xy_0_0, g_yy_0_0_0_xy_xz_0_0, g_yy_0_0_0_xy_yy_0_0, g_yy_0_0_0_xy_yz_0_0, g_yy_0_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xx_0_0[i] = -6.0 * g_xy_xx_0_0[i] * a_exp + 4.0 * g_xyyy_xx_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_0_0[i] = -6.0 * g_xy_xy_0_0[i] * a_exp + 4.0 * g_xyyy_xy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_0_0[i] = -6.0 * g_xy_xz_0_0[i] * a_exp + 4.0 * g_xyyy_xz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_0_0[i] = -6.0 * g_xy_yy_0_0[i] * a_exp + 4.0 * g_xyyy_yy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_0_0[i] = -6.0 * g_xy_yz_0_0[i] * a_exp + 4.0 * g_xyyy_yz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_0_0[i] = -6.0 * g_xy_zz_0_0[i] * a_exp + 4.0 * g_xyyy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xyyz_xx_0_0, g_xyyz_xy_0_0, g_xyyz_xz_0_0, g_xyyz_yy_0_0, g_xyyz_yz_0_0, g_xyyz_zz_0_0, g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_xz_zz_0_0, g_yy_0_0_0_xz_xx_0_0, g_yy_0_0_0_xz_xy_0_0, g_yy_0_0_0_xz_xz_0_0, g_yy_0_0_0_xz_yy_0_0, g_yy_0_0_0_xz_yz_0_0, g_yy_0_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xx_0_0[i] = -2.0 * g_xz_xx_0_0[i] * a_exp + 4.0 * g_xyyz_xx_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_0_0[i] = -2.0 * g_xz_xy_0_0[i] * a_exp + 4.0 * g_xyyz_xy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_0_0[i] = -2.0 * g_xz_xz_0_0[i] * a_exp + 4.0 * g_xyyz_xz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_0_0[i] = -2.0 * g_xz_yy_0_0[i] * a_exp + 4.0 * g_xyyz_yy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_0_0[i] = -2.0 * g_xz_yz_0_0[i] * a_exp + 4.0 * g_xyyz_yz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_0_0[i] = -2.0 * g_xz_zz_0_0[i] * a_exp + 4.0 * g_xyyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_yy_0_0_0_yy_xx_0_0, g_yy_0_0_0_yy_xy_0_0, g_yy_0_0_0_yy_xz_0_0, g_yy_0_0_0_yy_yy_0_0, g_yy_0_0_0_yy_yz_0_0, g_yy_0_0_0_yy_zz_0_0, g_yy_xx_0_0, g_yy_xy_0_0, g_yy_xz_0_0, g_yy_yy_0_0, g_yy_yz_0_0, g_yy_zz_0_0, g_yyyy_xx_0_0, g_yyyy_xy_0_0, g_yyyy_xz_0_0, g_yyyy_yy_0_0, g_yyyy_yz_0_0, g_yyyy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xx_0_0[i] = 2.0 * g_0_xx_0_0[i] - 10.0 * g_yy_xx_0_0[i] * a_exp + 4.0 * g_yyyy_xx_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_0_0[i] = 2.0 * g_0_xy_0_0[i] - 10.0 * g_yy_xy_0_0[i] * a_exp + 4.0 * g_yyyy_xy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_0_0[i] = 2.0 * g_0_xz_0_0[i] - 10.0 * g_yy_xz_0_0[i] * a_exp + 4.0 * g_yyyy_xz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_0_0[i] = 2.0 * g_0_yy_0_0[i] - 10.0 * g_yy_yy_0_0[i] * a_exp + 4.0 * g_yyyy_yy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_0_0[i] = 2.0 * g_0_yz_0_0[i] - 10.0 * g_yy_yz_0_0[i] * a_exp + 4.0 * g_yyyy_yz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_0_0[i] = 2.0 * g_0_zz_0_0[i] - 10.0 * g_yy_zz_0_0[i] * a_exp + 4.0 * g_yyyy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xx_0_0, g_yy_0_0_0_yz_xy_0_0, g_yy_0_0_0_yz_xz_0_0, g_yy_0_0_0_yz_yy_0_0, g_yy_0_0_0_yz_yz_0_0, g_yy_0_0_0_yz_zz_0_0, g_yyyz_xx_0_0, g_yyyz_xy_0_0, g_yyyz_xz_0_0, g_yyyz_yy_0_0, g_yyyz_yz_0_0, g_yyyz_zz_0_0, g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xx_0_0[i] = -6.0 * g_yz_xx_0_0[i] * a_exp + 4.0 * g_yyyz_xx_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_0_0[i] = -6.0 * g_yz_xy_0_0[i] * a_exp + 4.0 * g_yyyz_xy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_0_0[i] = -6.0 * g_yz_xz_0_0[i] * a_exp + 4.0 * g_yyyz_xz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_0_0[i] = -6.0 * g_yz_yy_0_0[i] * a_exp + 4.0 * g_yyyz_yy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_0_0[i] = -6.0 * g_yz_yz_0_0[i] * a_exp + 4.0 * g_yyyz_yz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_0_0[i] = -6.0 * g_yz_zz_0_0[i] * a_exp + 4.0 * g_yyyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xx_0_0, g_yy_0_0_0_zz_xy_0_0, g_yy_0_0_0_zz_xz_0_0, g_yy_0_0_0_zz_yy_0_0, g_yy_0_0_0_zz_yz_0_0, g_yy_0_0_0_zz_zz_0_0, g_yyzz_xx_0_0, g_yyzz_xy_0_0, g_yyzz_xz_0_0, g_yyzz_yy_0_0, g_yyzz_yz_0_0, g_yyzz_zz_0_0, g_zz_xx_0_0, g_zz_xy_0_0, g_zz_xz_0_0, g_zz_yy_0_0, g_zz_yz_0_0, g_zz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xx_0_0[i] = -2.0 * g_zz_xx_0_0[i] * a_exp + 4.0 * g_yyzz_xx_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_0_0[i] = -2.0 * g_zz_xy_0_0[i] * a_exp + 4.0 * g_yyzz_xy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_0_0[i] = -2.0 * g_zz_xz_0_0[i] * a_exp + 4.0 * g_yyzz_xz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_0_0[i] = -2.0 * g_zz_yy_0_0[i] * a_exp + 4.0 * g_yyzz_yy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_0_0[i] = -2.0 * g_zz_yz_0_0[i] * a_exp + 4.0 * g_yyzz_yz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_0_0[i] = -2.0 * g_zz_zz_0_0[i] * a_exp + 4.0 * g_yyzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xxyz_xx_0_0, g_xxyz_xy_0_0, g_xxyz_xz_0_0, g_xxyz_yy_0_0, g_xxyz_yz_0_0, g_xxyz_zz_0_0, g_yz_0_0_0_xx_xx_0_0, g_yz_0_0_0_xx_xy_0_0, g_yz_0_0_0_xx_xz_0_0, g_yz_0_0_0_xx_yy_0_0, g_yz_0_0_0_xx_yz_0_0, g_yz_0_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xx_0_0[i] = 4.0 * g_xxyz_xx_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_0_0[i] = 4.0 * g_xxyz_xy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_0_0[i] = 4.0 * g_xxyz_xz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_0_0[i] = 4.0 * g_xxyz_yy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_0_0[i] = 4.0 * g_xxyz_yz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_0_0[i] = 4.0 * g_xxyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_xyyz_xx_0_0, g_xyyz_xy_0_0, g_xyyz_xz_0_0, g_xyyz_yy_0_0, g_xyyz_yz_0_0, g_xyyz_zz_0_0, g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_xz_zz_0_0, g_yz_0_0_0_xy_xx_0_0, g_yz_0_0_0_xy_xy_0_0, g_yz_0_0_0_xy_xz_0_0, g_yz_0_0_0_xy_yy_0_0, g_yz_0_0_0_xy_yz_0_0, g_yz_0_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xx_0_0[i] = -2.0 * g_xz_xx_0_0[i] * a_exp + 4.0 * g_xyyz_xx_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_0_0[i] = -2.0 * g_xz_xy_0_0[i] * a_exp + 4.0 * g_xyyz_xy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_0_0[i] = -2.0 * g_xz_xz_0_0[i] * a_exp + 4.0 * g_xyyz_xz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_0_0[i] = -2.0 * g_xz_yy_0_0[i] * a_exp + 4.0 * g_xyyz_yy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_0_0[i] = -2.0 * g_xz_yz_0_0[i] * a_exp + 4.0 * g_xyyz_yz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_0_0[i] = -2.0 * g_xz_zz_0_0[i] * a_exp + 4.0 * g_xyyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_xy_zz_0_0, g_xyzz_xx_0_0, g_xyzz_xy_0_0, g_xyzz_xz_0_0, g_xyzz_yy_0_0, g_xyzz_yz_0_0, g_xyzz_zz_0_0, g_yz_0_0_0_xz_xx_0_0, g_yz_0_0_0_xz_xy_0_0, g_yz_0_0_0_xz_xz_0_0, g_yz_0_0_0_xz_yy_0_0, g_yz_0_0_0_xz_yz_0_0, g_yz_0_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xx_0_0[i] = -2.0 * g_xy_xx_0_0[i] * a_exp + 4.0 * g_xyzz_xx_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_0_0[i] = -2.0 * g_xy_xy_0_0[i] * a_exp + 4.0 * g_xyzz_xy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_0_0[i] = -2.0 * g_xy_xz_0_0[i] * a_exp + 4.0 * g_xyzz_xz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_0_0[i] = -2.0 * g_xy_yy_0_0[i] * a_exp + 4.0 * g_xyzz_yy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_0_0[i] = -2.0 * g_xy_yz_0_0[i] * a_exp + 4.0 * g_xyzz_yz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_0_0[i] = -2.0 * g_xy_zz_0_0[i] * a_exp + 4.0 * g_xyzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_yyyz_xx_0_0, g_yyyz_xy_0_0, g_yyyz_xz_0_0, g_yyyz_yy_0_0, g_yyyz_yz_0_0, g_yyyz_zz_0_0, g_yz_0_0_0_yy_xx_0_0, g_yz_0_0_0_yy_xy_0_0, g_yz_0_0_0_yy_xz_0_0, g_yz_0_0_0_yy_yy_0_0, g_yz_0_0_0_yy_yz_0_0, g_yz_0_0_0_yy_zz_0_0, g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xx_0_0[i] = -4.0 * g_yz_xx_0_0[i] * a_exp + 4.0 * g_yyyz_xx_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_0_0[i] = -4.0 * g_yz_xy_0_0[i] * a_exp + 4.0 * g_yyyz_xy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_0_0[i] = -4.0 * g_yz_xz_0_0[i] * a_exp + 4.0 * g_yyyz_xz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_0_0[i] = -4.0 * g_yz_yy_0_0[i] * a_exp + 4.0 * g_yyyz_yy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_0_0[i] = -4.0 * g_yz_yz_0_0[i] * a_exp + 4.0 * g_yyyz_yz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_0_0[i] = -4.0 * g_yz_zz_0_0[i] * a_exp + 4.0 * g_yyyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_yy_xx_0_0, g_yy_xy_0_0, g_yy_xz_0_0, g_yy_yy_0_0, g_yy_yz_0_0, g_yy_zz_0_0, g_yyzz_xx_0_0, g_yyzz_xy_0_0, g_yyzz_xz_0_0, g_yyzz_yy_0_0, g_yyzz_yz_0_0, g_yyzz_zz_0_0, g_yz_0_0_0_yz_xx_0_0, g_yz_0_0_0_yz_xy_0_0, g_yz_0_0_0_yz_xz_0_0, g_yz_0_0_0_yz_yy_0_0, g_yz_0_0_0_yz_yz_0_0, g_yz_0_0_0_yz_zz_0_0, g_zz_xx_0_0, g_zz_xy_0_0, g_zz_xz_0_0, g_zz_yy_0_0, g_zz_yz_0_0, g_zz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xx_0_0[i] = g_0_xx_0_0[i] - 2.0 * g_zz_xx_0_0[i] * a_exp - 2.0 * g_yy_xx_0_0[i] * a_exp + 4.0 * g_yyzz_xx_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_0_0[i] = g_0_xy_0_0[i] - 2.0 * g_zz_xy_0_0[i] * a_exp - 2.0 * g_yy_xy_0_0[i] * a_exp + 4.0 * g_yyzz_xy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_0_0[i] = g_0_xz_0_0[i] - 2.0 * g_zz_xz_0_0[i] * a_exp - 2.0 * g_yy_xz_0_0[i] * a_exp + 4.0 * g_yyzz_xz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_0_0[i] = g_0_yy_0_0[i] - 2.0 * g_zz_yy_0_0[i] * a_exp - 2.0 * g_yy_yy_0_0[i] * a_exp + 4.0 * g_yyzz_yy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_0_0[i] = g_0_yz_0_0[i] - 2.0 * g_zz_yz_0_0[i] * a_exp - 2.0 * g_yy_yz_0_0[i] * a_exp + 4.0 * g_yyzz_yz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_0_0[i] = g_0_zz_0_0[i] - 2.0 * g_zz_zz_0_0[i] * a_exp - 2.0 * g_yy_zz_0_0[i] * a_exp + 4.0 * g_yyzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xx_0_0, g_yz_0_0_0_zz_xy_0_0, g_yz_0_0_0_zz_xz_0_0, g_yz_0_0_0_zz_yy_0_0, g_yz_0_0_0_zz_yz_0_0, g_yz_0_0_0_zz_zz_0_0, g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_yz_zz_0_0, g_yzzz_xx_0_0, g_yzzz_xy_0_0, g_yzzz_xz_0_0, g_yzzz_yy_0_0, g_yzzz_yz_0_0, g_yzzz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xx_0_0[i] = -4.0 * g_yz_xx_0_0[i] * a_exp + 4.0 * g_yzzz_xx_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_0_0[i] = -4.0 * g_yz_xy_0_0[i] * a_exp + 4.0 * g_yzzz_xy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_0_0[i] = -4.0 * g_yz_xz_0_0[i] * a_exp + 4.0 * g_yzzz_xz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_0_0[i] = -4.0 * g_yz_yy_0_0[i] * a_exp + 4.0 * g_yzzz_yy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_0_0[i] = -4.0 * g_yz_yz_0_0[i] * a_exp + 4.0 * g_yzzz_yz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_0_0[i] = -4.0 * g_yz_zz_0_0[i] * a_exp + 4.0 * g_yzzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xx_xx_0_0, g_xx_xy_0_0, g_xx_xz_0_0, g_xx_yy_0_0, g_xx_yz_0_0, g_xx_zz_0_0, g_xxzz_xx_0_0, g_xxzz_xy_0_0, g_xxzz_xz_0_0, g_xxzz_yy_0_0, g_xxzz_yz_0_0, g_xxzz_zz_0_0, g_zz_0_0_0_xx_xx_0_0, g_zz_0_0_0_xx_xy_0_0, g_zz_0_0_0_xx_xz_0_0, g_zz_0_0_0_xx_yy_0_0, g_zz_0_0_0_xx_yz_0_0, g_zz_0_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xx_0_0[i] = -2.0 * g_xx_xx_0_0[i] * a_exp + 4.0 * g_xxzz_xx_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_0_0[i] = -2.0 * g_xx_xy_0_0[i] * a_exp + 4.0 * g_xxzz_xy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_0_0[i] = -2.0 * g_xx_xz_0_0[i] * a_exp + 4.0 * g_xxzz_xz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_0_0[i] = -2.0 * g_xx_yy_0_0[i] * a_exp + 4.0 * g_xxzz_yy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_0_0[i] = -2.0 * g_xx_yz_0_0[i] * a_exp + 4.0 * g_xxzz_yz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_0_0[i] = -2.0 * g_xx_zz_0_0[i] * a_exp + 4.0 * g_xxzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_xy_zz_0_0, g_xyzz_xx_0_0, g_xyzz_xy_0_0, g_xyzz_xz_0_0, g_xyzz_yy_0_0, g_xyzz_yz_0_0, g_xyzz_zz_0_0, g_zz_0_0_0_xy_xx_0_0, g_zz_0_0_0_xy_xy_0_0, g_zz_0_0_0_xy_xz_0_0, g_zz_0_0_0_xy_yy_0_0, g_zz_0_0_0_xy_yz_0_0, g_zz_0_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xx_0_0[i] = -2.0 * g_xy_xx_0_0[i] * a_exp + 4.0 * g_xyzz_xx_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_0_0[i] = -2.0 * g_xy_xy_0_0[i] * a_exp + 4.0 * g_xyzz_xy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_0_0[i] = -2.0 * g_xy_xz_0_0[i] * a_exp + 4.0 * g_xyzz_xz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_0_0[i] = -2.0 * g_xy_yy_0_0[i] * a_exp + 4.0 * g_xyzz_yy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_0_0[i] = -2.0 * g_xy_yz_0_0[i] * a_exp + 4.0 * g_xyzz_yz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_0_0[i] = -2.0 * g_xy_zz_0_0[i] * a_exp + 4.0 * g_xyzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_xz_zz_0_0, g_xzzz_xx_0_0, g_xzzz_xy_0_0, g_xzzz_xz_0_0, g_xzzz_yy_0_0, g_xzzz_yz_0_0, g_xzzz_zz_0_0, g_zz_0_0_0_xz_xx_0_0, g_zz_0_0_0_xz_xy_0_0, g_zz_0_0_0_xz_xz_0_0, g_zz_0_0_0_xz_yy_0_0, g_zz_0_0_0_xz_yz_0_0, g_zz_0_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xx_0_0[i] = -6.0 * g_xz_xx_0_0[i] * a_exp + 4.0 * g_xzzz_xx_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_0_0[i] = -6.0 * g_xz_xy_0_0[i] * a_exp + 4.0 * g_xzzz_xy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_0_0[i] = -6.0 * g_xz_xz_0_0[i] * a_exp + 4.0 * g_xzzz_xz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_0_0[i] = -6.0 * g_xz_yy_0_0[i] * a_exp + 4.0 * g_xzzz_yy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_0_0[i] = -6.0 * g_xz_yz_0_0[i] * a_exp + 4.0 * g_xzzz_yz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_0_0[i] = -6.0 * g_xz_zz_0_0[i] * a_exp + 4.0 * g_xzzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_yy_xx_0_0, g_yy_xy_0_0, g_yy_xz_0_0, g_yy_yy_0_0, g_yy_yz_0_0, g_yy_zz_0_0, g_yyzz_xx_0_0, g_yyzz_xy_0_0, g_yyzz_xz_0_0, g_yyzz_yy_0_0, g_yyzz_yz_0_0, g_yyzz_zz_0_0, g_zz_0_0_0_yy_xx_0_0, g_zz_0_0_0_yy_xy_0_0, g_zz_0_0_0_yy_xz_0_0, g_zz_0_0_0_yy_yy_0_0, g_zz_0_0_0_yy_yz_0_0, g_zz_0_0_0_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xx_0_0[i] = -2.0 * g_yy_xx_0_0[i] * a_exp + 4.0 * g_yyzz_xx_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_0_0[i] = -2.0 * g_yy_xy_0_0[i] * a_exp + 4.0 * g_yyzz_xy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_0_0[i] = -2.0 * g_yy_xz_0_0[i] * a_exp + 4.0 * g_yyzz_xz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_0_0[i] = -2.0 * g_yy_yy_0_0[i] * a_exp + 4.0 * g_yyzz_yy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_0_0[i] = -2.0 * g_yy_yz_0_0[i] * a_exp + 4.0 * g_yyzz_yz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_0_0[i] = -2.0 * g_yy_zz_0_0[i] * a_exp + 4.0 * g_yyzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_yz_zz_0_0, g_yzzz_xx_0_0, g_yzzz_xy_0_0, g_yzzz_xz_0_0, g_yzzz_yy_0_0, g_yzzz_yz_0_0, g_yzzz_zz_0_0, g_zz_0_0_0_yz_xx_0_0, g_zz_0_0_0_yz_xy_0_0, g_zz_0_0_0_yz_xz_0_0, g_zz_0_0_0_yz_yy_0_0, g_zz_0_0_0_yz_yz_0_0, g_zz_0_0_0_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xx_0_0[i] = -6.0 * g_yz_xx_0_0[i] * a_exp + 4.0 * g_yzzz_xx_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_0_0[i] = -6.0 * g_yz_xy_0_0[i] * a_exp + 4.0 * g_yzzz_xy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_0_0[i] = -6.0 * g_yz_xz_0_0[i] * a_exp + 4.0 * g_yzzz_xz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_0_0[i] = -6.0 * g_yz_yy_0_0[i] * a_exp + 4.0 * g_yzzz_yy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_0_0[i] = -6.0 * g_yz_yz_0_0[i] * a_exp + 4.0 * g_yzzz_yz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_0_0[i] = -6.0 * g_yz_zz_0_0[i] * a_exp + 4.0 * g_yzzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_zz_0_0_0_zz_xx_0_0, g_zz_0_0_0_zz_xy_0_0, g_zz_0_0_0_zz_xz_0_0, g_zz_0_0_0_zz_yy_0_0, g_zz_0_0_0_zz_yz_0_0, g_zz_0_0_0_zz_zz_0_0, g_zz_xx_0_0, g_zz_xy_0_0, g_zz_xz_0_0, g_zz_yy_0_0, g_zz_yz_0_0, g_zz_zz_0_0, g_zzzz_xx_0_0, g_zzzz_xy_0_0, g_zzzz_xz_0_0, g_zzzz_yy_0_0, g_zzzz_yz_0_0, g_zzzz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xx_0_0[i] = 2.0 * g_0_xx_0_0[i] - 10.0 * g_zz_xx_0_0[i] * a_exp + 4.0 * g_zzzz_xx_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_0_0[i] = 2.0 * g_0_xy_0_0[i] - 10.0 * g_zz_xy_0_0[i] * a_exp + 4.0 * g_zzzz_xy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_0_0[i] = 2.0 * g_0_xz_0_0[i] - 10.0 * g_zz_xz_0_0[i] * a_exp + 4.0 * g_zzzz_xz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_0_0[i] = 2.0 * g_0_yy_0_0[i] - 10.0 * g_zz_yy_0_0[i] * a_exp + 4.0 * g_zzzz_yy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_0_0[i] = 2.0 * g_0_yz_0_0[i] - 10.0 * g_zz_yz_0_0[i] * a_exp + 4.0 * g_zzzz_yz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_0_0[i] = 2.0 * g_0_zz_0_0[i] - 10.0 * g_zz_zz_0_0[i] * a_exp + 4.0 * g_zzzz_zz_0_0[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

