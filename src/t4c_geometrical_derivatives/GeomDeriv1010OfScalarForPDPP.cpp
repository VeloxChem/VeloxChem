#include "GeomDeriv1010OfScalarForPDPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_pdpp_0(CSimdArray<double>& buffer_1010_pdpp,
                     const CSimdArray<double>& buffer_sdsp,
                     const CSimdArray<double>& buffer_sddp,
                     const CSimdArray<double>& buffer_ddsp,
                     const CSimdArray<double>& buffer_dddp,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_pdpp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sdsp

    auto g_0_xx_0_x = buffer_sdsp[0];

    auto g_0_xx_0_y = buffer_sdsp[1];

    auto g_0_xx_0_z = buffer_sdsp[2];

    auto g_0_xy_0_x = buffer_sdsp[3];

    auto g_0_xy_0_y = buffer_sdsp[4];

    auto g_0_xy_0_z = buffer_sdsp[5];

    auto g_0_xz_0_x = buffer_sdsp[6];

    auto g_0_xz_0_y = buffer_sdsp[7];

    auto g_0_xz_0_z = buffer_sdsp[8];

    auto g_0_yy_0_x = buffer_sdsp[9];

    auto g_0_yy_0_y = buffer_sdsp[10];

    auto g_0_yy_0_z = buffer_sdsp[11];

    auto g_0_yz_0_x = buffer_sdsp[12];

    auto g_0_yz_0_y = buffer_sdsp[13];

    auto g_0_yz_0_z = buffer_sdsp[14];

    auto g_0_zz_0_x = buffer_sdsp[15];

    auto g_0_zz_0_y = buffer_sdsp[16];

    auto g_0_zz_0_z = buffer_sdsp[17];

    /// Set up components of auxilary buffer : buffer_sddp

    auto g_0_xx_xx_x = buffer_sddp[0];

    auto g_0_xx_xx_y = buffer_sddp[1];

    auto g_0_xx_xx_z = buffer_sddp[2];

    auto g_0_xx_xy_x = buffer_sddp[3];

    auto g_0_xx_xy_y = buffer_sddp[4];

    auto g_0_xx_xy_z = buffer_sddp[5];

    auto g_0_xx_xz_x = buffer_sddp[6];

    auto g_0_xx_xz_y = buffer_sddp[7];

    auto g_0_xx_xz_z = buffer_sddp[8];

    auto g_0_xx_yy_x = buffer_sddp[9];

    auto g_0_xx_yy_y = buffer_sddp[10];

    auto g_0_xx_yy_z = buffer_sddp[11];

    auto g_0_xx_yz_x = buffer_sddp[12];

    auto g_0_xx_yz_y = buffer_sddp[13];

    auto g_0_xx_yz_z = buffer_sddp[14];

    auto g_0_xx_zz_x = buffer_sddp[15];

    auto g_0_xx_zz_y = buffer_sddp[16];

    auto g_0_xx_zz_z = buffer_sddp[17];

    auto g_0_xy_xx_x = buffer_sddp[18];

    auto g_0_xy_xx_y = buffer_sddp[19];

    auto g_0_xy_xx_z = buffer_sddp[20];

    auto g_0_xy_xy_x = buffer_sddp[21];

    auto g_0_xy_xy_y = buffer_sddp[22];

    auto g_0_xy_xy_z = buffer_sddp[23];

    auto g_0_xy_xz_x = buffer_sddp[24];

    auto g_0_xy_xz_y = buffer_sddp[25];

    auto g_0_xy_xz_z = buffer_sddp[26];

    auto g_0_xy_yy_x = buffer_sddp[27];

    auto g_0_xy_yy_y = buffer_sddp[28];

    auto g_0_xy_yy_z = buffer_sddp[29];

    auto g_0_xy_yz_x = buffer_sddp[30];

    auto g_0_xy_yz_y = buffer_sddp[31];

    auto g_0_xy_yz_z = buffer_sddp[32];

    auto g_0_xy_zz_x = buffer_sddp[33];

    auto g_0_xy_zz_y = buffer_sddp[34];

    auto g_0_xy_zz_z = buffer_sddp[35];

    auto g_0_xz_xx_x = buffer_sddp[36];

    auto g_0_xz_xx_y = buffer_sddp[37];

    auto g_0_xz_xx_z = buffer_sddp[38];

    auto g_0_xz_xy_x = buffer_sddp[39];

    auto g_0_xz_xy_y = buffer_sddp[40];

    auto g_0_xz_xy_z = buffer_sddp[41];

    auto g_0_xz_xz_x = buffer_sddp[42];

    auto g_0_xz_xz_y = buffer_sddp[43];

    auto g_0_xz_xz_z = buffer_sddp[44];

    auto g_0_xz_yy_x = buffer_sddp[45];

    auto g_0_xz_yy_y = buffer_sddp[46];

    auto g_0_xz_yy_z = buffer_sddp[47];

    auto g_0_xz_yz_x = buffer_sddp[48];

    auto g_0_xz_yz_y = buffer_sddp[49];

    auto g_0_xz_yz_z = buffer_sddp[50];

    auto g_0_xz_zz_x = buffer_sddp[51];

    auto g_0_xz_zz_y = buffer_sddp[52];

    auto g_0_xz_zz_z = buffer_sddp[53];

    auto g_0_yy_xx_x = buffer_sddp[54];

    auto g_0_yy_xx_y = buffer_sddp[55];

    auto g_0_yy_xx_z = buffer_sddp[56];

    auto g_0_yy_xy_x = buffer_sddp[57];

    auto g_0_yy_xy_y = buffer_sddp[58];

    auto g_0_yy_xy_z = buffer_sddp[59];

    auto g_0_yy_xz_x = buffer_sddp[60];

    auto g_0_yy_xz_y = buffer_sddp[61];

    auto g_0_yy_xz_z = buffer_sddp[62];

    auto g_0_yy_yy_x = buffer_sddp[63];

    auto g_0_yy_yy_y = buffer_sddp[64];

    auto g_0_yy_yy_z = buffer_sddp[65];

    auto g_0_yy_yz_x = buffer_sddp[66];

    auto g_0_yy_yz_y = buffer_sddp[67];

    auto g_0_yy_yz_z = buffer_sddp[68];

    auto g_0_yy_zz_x = buffer_sddp[69];

    auto g_0_yy_zz_y = buffer_sddp[70];

    auto g_0_yy_zz_z = buffer_sddp[71];

    auto g_0_yz_xx_x = buffer_sddp[72];

    auto g_0_yz_xx_y = buffer_sddp[73];

    auto g_0_yz_xx_z = buffer_sddp[74];

    auto g_0_yz_xy_x = buffer_sddp[75];

    auto g_0_yz_xy_y = buffer_sddp[76];

    auto g_0_yz_xy_z = buffer_sddp[77];

    auto g_0_yz_xz_x = buffer_sddp[78];

    auto g_0_yz_xz_y = buffer_sddp[79];

    auto g_0_yz_xz_z = buffer_sddp[80];

    auto g_0_yz_yy_x = buffer_sddp[81];

    auto g_0_yz_yy_y = buffer_sddp[82];

    auto g_0_yz_yy_z = buffer_sddp[83];

    auto g_0_yz_yz_x = buffer_sddp[84];

    auto g_0_yz_yz_y = buffer_sddp[85];

    auto g_0_yz_yz_z = buffer_sddp[86];

    auto g_0_yz_zz_x = buffer_sddp[87];

    auto g_0_yz_zz_y = buffer_sddp[88];

    auto g_0_yz_zz_z = buffer_sddp[89];

    auto g_0_zz_xx_x = buffer_sddp[90];

    auto g_0_zz_xx_y = buffer_sddp[91];

    auto g_0_zz_xx_z = buffer_sddp[92];

    auto g_0_zz_xy_x = buffer_sddp[93];

    auto g_0_zz_xy_y = buffer_sddp[94];

    auto g_0_zz_xy_z = buffer_sddp[95];

    auto g_0_zz_xz_x = buffer_sddp[96];

    auto g_0_zz_xz_y = buffer_sddp[97];

    auto g_0_zz_xz_z = buffer_sddp[98];

    auto g_0_zz_yy_x = buffer_sddp[99];

    auto g_0_zz_yy_y = buffer_sddp[100];

    auto g_0_zz_yy_z = buffer_sddp[101];

    auto g_0_zz_yz_x = buffer_sddp[102];

    auto g_0_zz_yz_y = buffer_sddp[103];

    auto g_0_zz_yz_z = buffer_sddp[104];

    auto g_0_zz_zz_x = buffer_sddp[105];

    auto g_0_zz_zz_y = buffer_sddp[106];

    auto g_0_zz_zz_z = buffer_sddp[107];

    /// Set up components of auxilary buffer : buffer_ddsp

    auto g_xx_xx_0_x = buffer_ddsp[0];

    auto g_xx_xx_0_y = buffer_ddsp[1];

    auto g_xx_xx_0_z = buffer_ddsp[2];

    auto g_xx_xy_0_x = buffer_ddsp[3];

    auto g_xx_xy_0_y = buffer_ddsp[4];

    auto g_xx_xy_0_z = buffer_ddsp[5];

    auto g_xx_xz_0_x = buffer_ddsp[6];

    auto g_xx_xz_0_y = buffer_ddsp[7];

    auto g_xx_xz_0_z = buffer_ddsp[8];

    auto g_xx_yy_0_x = buffer_ddsp[9];

    auto g_xx_yy_0_y = buffer_ddsp[10];

    auto g_xx_yy_0_z = buffer_ddsp[11];

    auto g_xx_yz_0_x = buffer_ddsp[12];

    auto g_xx_yz_0_y = buffer_ddsp[13];

    auto g_xx_yz_0_z = buffer_ddsp[14];

    auto g_xx_zz_0_x = buffer_ddsp[15];

    auto g_xx_zz_0_y = buffer_ddsp[16];

    auto g_xx_zz_0_z = buffer_ddsp[17];

    auto g_xy_xx_0_x = buffer_ddsp[18];

    auto g_xy_xx_0_y = buffer_ddsp[19];

    auto g_xy_xx_0_z = buffer_ddsp[20];

    auto g_xy_xy_0_x = buffer_ddsp[21];

    auto g_xy_xy_0_y = buffer_ddsp[22];

    auto g_xy_xy_0_z = buffer_ddsp[23];

    auto g_xy_xz_0_x = buffer_ddsp[24];

    auto g_xy_xz_0_y = buffer_ddsp[25];

    auto g_xy_xz_0_z = buffer_ddsp[26];

    auto g_xy_yy_0_x = buffer_ddsp[27];

    auto g_xy_yy_0_y = buffer_ddsp[28];

    auto g_xy_yy_0_z = buffer_ddsp[29];

    auto g_xy_yz_0_x = buffer_ddsp[30];

    auto g_xy_yz_0_y = buffer_ddsp[31];

    auto g_xy_yz_0_z = buffer_ddsp[32];

    auto g_xy_zz_0_x = buffer_ddsp[33];

    auto g_xy_zz_0_y = buffer_ddsp[34];

    auto g_xy_zz_0_z = buffer_ddsp[35];

    auto g_xz_xx_0_x = buffer_ddsp[36];

    auto g_xz_xx_0_y = buffer_ddsp[37];

    auto g_xz_xx_0_z = buffer_ddsp[38];

    auto g_xz_xy_0_x = buffer_ddsp[39];

    auto g_xz_xy_0_y = buffer_ddsp[40];

    auto g_xz_xy_0_z = buffer_ddsp[41];

    auto g_xz_xz_0_x = buffer_ddsp[42];

    auto g_xz_xz_0_y = buffer_ddsp[43];

    auto g_xz_xz_0_z = buffer_ddsp[44];

    auto g_xz_yy_0_x = buffer_ddsp[45];

    auto g_xz_yy_0_y = buffer_ddsp[46];

    auto g_xz_yy_0_z = buffer_ddsp[47];

    auto g_xz_yz_0_x = buffer_ddsp[48];

    auto g_xz_yz_0_y = buffer_ddsp[49];

    auto g_xz_yz_0_z = buffer_ddsp[50];

    auto g_xz_zz_0_x = buffer_ddsp[51];

    auto g_xz_zz_0_y = buffer_ddsp[52];

    auto g_xz_zz_0_z = buffer_ddsp[53];

    auto g_yy_xx_0_x = buffer_ddsp[54];

    auto g_yy_xx_0_y = buffer_ddsp[55];

    auto g_yy_xx_0_z = buffer_ddsp[56];

    auto g_yy_xy_0_x = buffer_ddsp[57];

    auto g_yy_xy_0_y = buffer_ddsp[58];

    auto g_yy_xy_0_z = buffer_ddsp[59];

    auto g_yy_xz_0_x = buffer_ddsp[60];

    auto g_yy_xz_0_y = buffer_ddsp[61];

    auto g_yy_xz_0_z = buffer_ddsp[62];

    auto g_yy_yy_0_x = buffer_ddsp[63];

    auto g_yy_yy_0_y = buffer_ddsp[64];

    auto g_yy_yy_0_z = buffer_ddsp[65];

    auto g_yy_yz_0_x = buffer_ddsp[66];

    auto g_yy_yz_0_y = buffer_ddsp[67];

    auto g_yy_yz_0_z = buffer_ddsp[68];

    auto g_yy_zz_0_x = buffer_ddsp[69];

    auto g_yy_zz_0_y = buffer_ddsp[70];

    auto g_yy_zz_0_z = buffer_ddsp[71];

    auto g_yz_xx_0_x = buffer_ddsp[72];

    auto g_yz_xx_0_y = buffer_ddsp[73];

    auto g_yz_xx_0_z = buffer_ddsp[74];

    auto g_yz_xy_0_x = buffer_ddsp[75];

    auto g_yz_xy_0_y = buffer_ddsp[76];

    auto g_yz_xy_0_z = buffer_ddsp[77];

    auto g_yz_xz_0_x = buffer_ddsp[78];

    auto g_yz_xz_0_y = buffer_ddsp[79];

    auto g_yz_xz_0_z = buffer_ddsp[80];

    auto g_yz_yy_0_x = buffer_ddsp[81];

    auto g_yz_yy_0_y = buffer_ddsp[82];

    auto g_yz_yy_0_z = buffer_ddsp[83];

    auto g_yz_yz_0_x = buffer_ddsp[84];

    auto g_yz_yz_0_y = buffer_ddsp[85];

    auto g_yz_yz_0_z = buffer_ddsp[86];

    auto g_yz_zz_0_x = buffer_ddsp[87];

    auto g_yz_zz_0_y = buffer_ddsp[88];

    auto g_yz_zz_0_z = buffer_ddsp[89];

    auto g_zz_xx_0_x = buffer_ddsp[90];

    auto g_zz_xx_0_y = buffer_ddsp[91];

    auto g_zz_xx_0_z = buffer_ddsp[92];

    auto g_zz_xy_0_x = buffer_ddsp[93];

    auto g_zz_xy_0_y = buffer_ddsp[94];

    auto g_zz_xy_0_z = buffer_ddsp[95];

    auto g_zz_xz_0_x = buffer_ddsp[96];

    auto g_zz_xz_0_y = buffer_ddsp[97];

    auto g_zz_xz_0_z = buffer_ddsp[98];

    auto g_zz_yy_0_x = buffer_ddsp[99];

    auto g_zz_yy_0_y = buffer_ddsp[100];

    auto g_zz_yy_0_z = buffer_ddsp[101];

    auto g_zz_yz_0_x = buffer_ddsp[102];

    auto g_zz_yz_0_y = buffer_ddsp[103];

    auto g_zz_yz_0_z = buffer_ddsp[104];

    auto g_zz_zz_0_x = buffer_ddsp[105];

    auto g_zz_zz_0_y = buffer_ddsp[106];

    auto g_zz_zz_0_z = buffer_ddsp[107];

    /// Set up components of auxilary buffer : buffer_dddp

    auto g_xx_xx_xx_x = buffer_dddp[0];

    auto g_xx_xx_xx_y = buffer_dddp[1];

    auto g_xx_xx_xx_z = buffer_dddp[2];

    auto g_xx_xx_xy_x = buffer_dddp[3];

    auto g_xx_xx_xy_y = buffer_dddp[4];

    auto g_xx_xx_xy_z = buffer_dddp[5];

    auto g_xx_xx_xz_x = buffer_dddp[6];

    auto g_xx_xx_xz_y = buffer_dddp[7];

    auto g_xx_xx_xz_z = buffer_dddp[8];

    auto g_xx_xx_yy_x = buffer_dddp[9];

    auto g_xx_xx_yy_y = buffer_dddp[10];

    auto g_xx_xx_yy_z = buffer_dddp[11];

    auto g_xx_xx_yz_x = buffer_dddp[12];

    auto g_xx_xx_yz_y = buffer_dddp[13];

    auto g_xx_xx_yz_z = buffer_dddp[14];

    auto g_xx_xx_zz_x = buffer_dddp[15];

    auto g_xx_xx_zz_y = buffer_dddp[16];

    auto g_xx_xx_zz_z = buffer_dddp[17];

    auto g_xx_xy_xx_x = buffer_dddp[18];

    auto g_xx_xy_xx_y = buffer_dddp[19];

    auto g_xx_xy_xx_z = buffer_dddp[20];

    auto g_xx_xy_xy_x = buffer_dddp[21];

    auto g_xx_xy_xy_y = buffer_dddp[22];

    auto g_xx_xy_xy_z = buffer_dddp[23];

    auto g_xx_xy_xz_x = buffer_dddp[24];

    auto g_xx_xy_xz_y = buffer_dddp[25];

    auto g_xx_xy_xz_z = buffer_dddp[26];

    auto g_xx_xy_yy_x = buffer_dddp[27];

    auto g_xx_xy_yy_y = buffer_dddp[28];

    auto g_xx_xy_yy_z = buffer_dddp[29];

    auto g_xx_xy_yz_x = buffer_dddp[30];

    auto g_xx_xy_yz_y = buffer_dddp[31];

    auto g_xx_xy_yz_z = buffer_dddp[32];

    auto g_xx_xy_zz_x = buffer_dddp[33];

    auto g_xx_xy_zz_y = buffer_dddp[34];

    auto g_xx_xy_zz_z = buffer_dddp[35];

    auto g_xx_xz_xx_x = buffer_dddp[36];

    auto g_xx_xz_xx_y = buffer_dddp[37];

    auto g_xx_xz_xx_z = buffer_dddp[38];

    auto g_xx_xz_xy_x = buffer_dddp[39];

    auto g_xx_xz_xy_y = buffer_dddp[40];

    auto g_xx_xz_xy_z = buffer_dddp[41];

    auto g_xx_xz_xz_x = buffer_dddp[42];

    auto g_xx_xz_xz_y = buffer_dddp[43];

    auto g_xx_xz_xz_z = buffer_dddp[44];

    auto g_xx_xz_yy_x = buffer_dddp[45];

    auto g_xx_xz_yy_y = buffer_dddp[46];

    auto g_xx_xz_yy_z = buffer_dddp[47];

    auto g_xx_xz_yz_x = buffer_dddp[48];

    auto g_xx_xz_yz_y = buffer_dddp[49];

    auto g_xx_xz_yz_z = buffer_dddp[50];

    auto g_xx_xz_zz_x = buffer_dddp[51];

    auto g_xx_xz_zz_y = buffer_dddp[52];

    auto g_xx_xz_zz_z = buffer_dddp[53];

    auto g_xx_yy_xx_x = buffer_dddp[54];

    auto g_xx_yy_xx_y = buffer_dddp[55];

    auto g_xx_yy_xx_z = buffer_dddp[56];

    auto g_xx_yy_xy_x = buffer_dddp[57];

    auto g_xx_yy_xy_y = buffer_dddp[58];

    auto g_xx_yy_xy_z = buffer_dddp[59];

    auto g_xx_yy_xz_x = buffer_dddp[60];

    auto g_xx_yy_xz_y = buffer_dddp[61];

    auto g_xx_yy_xz_z = buffer_dddp[62];

    auto g_xx_yy_yy_x = buffer_dddp[63];

    auto g_xx_yy_yy_y = buffer_dddp[64];

    auto g_xx_yy_yy_z = buffer_dddp[65];

    auto g_xx_yy_yz_x = buffer_dddp[66];

    auto g_xx_yy_yz_y = buffer_dddp[67];

    auto g_xx_yy_yz_z = buffer_dddp[68];

    auto g_xx_yy_zz_x = buffer_dddp[69];

    auto g_xx_yy_zz_y = buffer_dddp[70];

    auto g_xx_yy_zz_z = buffer_dddp[71];

    auto g_xx_yz_xx_x = buffer_dddp[72];

    auto g_xx_yz_xx_y = buffer_dddp[73];

    auto g_xx_yz_xx_z = buffer_dddp[74];

    auto g_xx_yz_xy_x = buffer_dddp[75];

    auto g_xx_yz_xy_y = buffer_dddp[76];

    auto g_xx_yz_xy_z = buffer_dddp[77];

    auto g_xx_yz_xz_x = buffer_dddp[78];

    auto g_xx_yz_xz_y = buffer_dddp[79];

    auto g_xx_yz_xz_z = buffer_dddp[80];

    auto g_xx_yz_yy_x = buffer_dddp[81];

    auto g_xx_yz_yy_y = buffer_dddp[82];

    auto g_xx_yz_yy_z = buffer_dddp[83];

    auto g_xx_yz_yz_x = buffer_dddp[84];

    auto g_xx_yz_yz_y = buffer_dddp[85];

    auto g_xx_yz_yz_z = buffer_dddp[86];

    auto g_xx_yz_zz_x = buffer_dddp[87];

    auto g_xx_yz_zz_y = buffer_dddp[88];

    auto g_xx_yz_zz_z = buffer_dddp[89];

    auto g_xx_zz_xx_x = buffer_dddp[90];

    auto g_xx_zz_xx_y = buffer_dddp[91];

    auto g_xx_zz_xx_z = buffer_dddp[92];

    auto g_xx_zz_xy_x = buffer_dddp[93];

    auto g_xx_zz_xy_y = buffer_dddp[94];

    auto g_xx_zz_xy_z = buffer_dddp[95];

    auto g_xx_zz_xz_x = buffer_dddp[96];

    auto g_xx_zz_xz_y = buffer_dddp[97];

    auto g_xx_zz_xz_z = buffer_dddp[98];

    auto g_xx_zz_yy_x = buffer_dddp[99];

    auto g_xx_zz_yy_y = buffer_dddp[100];

    auto g_xx_zz_yy_z = buffer_dddp[101];

    auto g_xx_zz_yz_x = buffer_dddp[102];

    auto g_xx_zz_yz_y = buffer_dddp[103];

    auto g_xx_zz_yz_z = buffer_dddp[104];

    auto g_xx_zz_zz_x = buffer_dddp[105];

    auto g_xx_zz_zz_y = buffer_dddp[106];

    auto g_xx_zz_zz_z = buffer_dddp[107];

    auto g_xy_xx_xx_x = buffer_dddp[108];

    auto g_xy_xx_xx_y = buffer_dddp[109];

    auto g_xy_xx_xx_z = buffer_dddp[110];

    auto g_xy_xx_xy_x = buffer_dddp[111];

    auto g_xy_xx_xy_y = buffer_dddp[112];

    auto g_xy_xx_xy_z = buffer_dddp[113];

    auto g_xy_xx_xz_x = buffer_dddp[114];

    auto g_xy_xx_xz_y = buffer_dddp[115];

    auto g_xy_xx_xz_z = buffer_dddp[116];

    auto g_xy_xx_yy_x = buffer_dddp[117];

    auto g_xy_xx_yy_y = buffer_dddp[118];

    auto g_xy_xx_yy_z = buffer_dddp[119];

    auto g_xy_xx_yz_x = buffer_dddp[120];

    auto g_xy_xx_yz_y = buffer_dddp[121];

    auto g_xy_xx_yz_z = buffer_dddp[122];

    auto g_xy_xx_zz_x = buffer_dddp[123];

    auto g_xy_xx_zz_y = buffer_dddp[124];

    auto g_xy_xx_zz_z = buffer_dddp[125];

    auto g_xy_xy_xx_x = buffer_dddp[126];

    auto g_xy_xy_xx_y = buffer_dddp[127];

    auto g_xy_xy_xx_z = buffer_dddp[128];

    auto g_xy_xy_xy_x = buffer_dddp[129];

    auto g_xy_xy_xy_y = buffer_dddp[130];

    auto g_xy_xy_xy_z = buffer_dddp[131];

    auto g_xy_xy_xz_x = buffer_dddp[132];

    auto g_xy_xy_xz_y = buffer_dddp[133];

    auto g_xy_xy_xz_z = buffer_dddp[134];

    auto g_xy_xy_yy_x = buffer_dddp[135];

    auto g_xy_xy_yy_y = buffer_dddp[136];

    auto g_xy_xy_yy_z = buffer_dddp[137];

    auto g_xy_xy_yz_x = buffer_dddp[138];

    auto g_xy_xy_yz_y = buffer_dddp[139];

    auto g_xy_xy_yz_z = buffer_dddp[140];

    auto g_xy_xy_zz_x = buffer_dddp[141];

    auto g_xy_xy_zz_y = buffer_dddp[142];

    auto g_xy_xy_zz_z = buffer_dddp[143];

    auto g_xy_xz_xx_x = buffer_dddp[144];

    auto g_xy_xz_xx_y = buffer_dddp[145];

    auto g_xy_xz_xx_z = buffer_dddp[146];

    auto g_xy_xz_xy_x = buffer_dddp[147];

    auto g_xy_xz_xy_y = buffer_dddp[148];

    auto g_xy_xz_xy_z = buffer_dddp[149];

    auto g_xy_xz_xz_x = buffer_dddp[150];

    auto g_xy_xz_xz_y = buffer_dddp[151];

    auto g_xy_xz_xz_z = buffer_dddp[152];

    auto g_xy_xz_yy_x = buffer_dddp[153];

    auto g_xy_xz_yy_y = buffer_dddp[154];

    auto g_xy_xz_yy_z = buffer_dddp[155];

    auto g_xy_xz_yz_x = buffer_dddp[156];

    auto g_xy_xz_yz_y = buffer_dddp[157];

    auto g_xy_xz_yz_z = buffer_dddp[158];

    auto g_xy_xz_zz_x = buffer_dddp[159];

    auto g_xy_xz_zz_y = buffer_dddp[160];

    auto g_xy_xz_zz_z = buffer_dddp[161];

    auto g_xy_yy_xx_x = buffer_dddp[162];

    auto g_xy_yy_xx_y = buffer_dddp[163];

    auto g_xy_yy_xx_z = buffer_dddp[164];

    auto g_xy_yy_xy_x = buffer_dddp[165];

    auto g_xy_yy_xy_y = buffer_dddp[166];

    auto g_xy_yy_xy_z = buffer_dddp[167];

    auto g_xy_yy_xz_x = buffer_dddp[168];

    auto g_xy_yy_xz_y = buffer_dddp[169];

    auto g_xy_yy_xz_z = buffer_dddp[170];

    auto g_xy_yy_yy_x = buffer_dddp[171];

    auto g_xy_yy_yy_y = buffer_dddp[172];

    auto g_xy_yy_yy_z = buffer_dddp[173];

    auto g_xy_yy_yz_x = buffer_dddp[174];

    auto g_xy_yy_yz_y = buffer_dddp[175];

    auto g_xy_yy_yz_z = buffer_dddp[176];

    auto g_xy_yy_zz_x = buffer_dddp[177];

    auto g_xy_yy_zz_y = buffer_dddp[178];

    auto g_xy_yy_zz_z = buffer_dddp[179];

    auto g_xy_yz_xx_x = buffer_dddp[180];

    auto g_xy_yz_xx_y = buffer_dddp[181];

    auto g_xy_yz_xx_z = buffer_dddp[182];

    auto g_xy_yz_xy_x = buffer_dddp[183];

    auto g_xy_yz_xy_y = buffer_dddp[184];

    auto g_xy_yz_xy_z = buffer_dddp[185];

    auto g_xy_yz_xz_x = buffer_dddp[186];

    auto g_xy_yz_xz_y = buffer_dddp[187];

    auto g_xy_yz_xz_z = buffer_dddp[188];

    auto g_xy_yz_yy_x = buffer_dddp[189];

    auto g_xy_yz_yy_y = buffer_dddp[190];

    auto g_xy_yz_yy_z = buffer_dddp[191];

    auto g_xy_yz_yz_x = buffer_dddp[192];

    auto g_xy_yz_yz_y = buffer_dddp[193];

    auto g_xy_yz_yz_z = buffer_dddp[194];

    auto g_xy_yz_zz_x = buffer_dddp[195];

    auto g_xy_yz_zz_y = buffer_dddp[196];

    auto g_xy_yz_zz_z = buffer_dddp[197];

    auto g_xy_zz_xx_x = buffer_dddp[198];

    auto g_xy_zz_xx_y = buffer_dddp[199];

    auto g_xy_zz_xx_z = buffer_dddp[200];

    auto g_xy_zz_xy_x = buffer_dddp[201];

    auto g_xy_zz_xy_y = buffer_dddp[202];

    auto g_xy_zz_xy_z = buffer_dddp[203];

    auto g_xy_zz_xz_x = buffer_dddp[204];

    auto g_xy_zz_xz_y = buffer_dddp[205];

    auto g_xy_zz_xz_z = buffer_dddp[206];

    auto g_xy_zz_yy_x = buffer_dddp[207];

    auto g_xy_zz_yy_y = buffer_dddp[208];

    auto g_xy_zz_yy_z = buffer_dddp[209];

    auto g_xy_zz_yz_x = buffer_dddp[210];

    auto g_xy_zz_yz_y = buffer_dddp[211];

    auto g_xy_zz_yz_z = buffer_dddp[212];

    auto g_xy_zz_zz_x = buffer_dddp[213];

    auto g_xy_zz_zz_y = buffer_dddp[214];

    auto g_xy_zz_zz_z = buffer_dddp[215];

    auto g_xz_xx_xx_x = buffer_dddp[216];

    auto g_xz_xx_xx_y = buffer_dddp[217];

    auto g_xz_xx_xx_z = buffer_dddp[218];

    auto g_xz_xx_xy_x = buffer_dddp[219];

    auto g_xz_xx_xy_y = buffer_dddp[220];

    auto g_xz_xx_xy_z = buffer_dddp[221];

    auto g_xz_xx_xz_x = buffer_dddp[222];

    auto g_xz_xx_xz_y = buffer_dddp[223];

    auto g_xz_xx_xz_z = buffer_dddp[224];

    auto g_xz_xx_yy_x = buffer_dddp[225];

    auto g_xz_xx_yy_y = buffer_dddp[226];

    auto g_xz_xx_yy_z = buffer_dddp[227];

    auto g_xz_xx_yz_x = buffer_dddp[228];

    auto g_xz_xx_yz_y = buffer_dddp[229];

    auto g_xz_xx_yz_z = buffer_dddp[230];

    auto g_xz_xx_zz_x = buffer_dddp[231];

    auto g_xz_xx_zz_y = buffer_dddp[232];

    auto g_xz_xx_zz_z = buffer_dddp[233];

    auto g_xz_xy_xx_x = buffer_dddp[234];

    auto g_xz_xy_xx_y = buffer_dddp[235];

    auto g_xz_xy_xx_z = buffer_dddp[236];

    auto g_xz_xy_xy_x = buffer_dddp[237];

    auto g_xz_xy_xy_y = buffer_dddp[238];

    auto g_xz_xy_xy_z = buffer_dddp[239];

    auto g_xz_xy_xz_x = buffer_dddp[240];

    auto g_xz_xy_xz_y = buffer_dddp[241];

    auto g_xz_xy_xz_z = buffer_dddp[242];

    auto g_xz_xy_yy_x = buffer_dddp[243];

    auto g_xz_xy_yy_y = buffer_dddp[244];

    auto g_xz_xy_yy_z = buffer_dddp[245];

    auto g_xz_xy_yz_x = buffer_dddp[246];

    auto g_xz_xy_yz_y = buffer_dddp[247];

    auto g_xz_xy_yz_z = buffer_dddp[248];

    auto g_xz_xy_zz_x = buffer_dddp[249];

    auto g_xz_xy_zz_y = buffer_dddp[250];

    auto g_xz_xy_zz_z = buffer_dddp[251];

    auto g_xz_xz_xx_x = buffer_dddp[252];

    auto g_xz_xz_xx_y = buffer_dddp[253];

    auto g_xz_xz_xx_z = buffer_dddp[254];

    auto g_xz_xz_xy_x = buffer_dddp[255];

    auto g_xz_xz_xy_y = buffer_dddp[256];

    auto g_xz_xz_xy_z = buffer_dddp[257];

    auto g_xz_xz_xz_x = buffer_dddp[258];

    auto g_xz_xz_xz_y = buffer_dddp[259];

    auto g_xz_xz_xz_z = buffer_dddp[260];

    auto g_xz_xz_yy_x = buffer_dddp[261];

    auto g_xz_xz_yy_y = buffer_dddp[262];

    auto g_xz_xz_yy_z = buffer_dddp[263];

    auto g_xz_xz_yz_x = buffer_dddp[264];

    auto g_xz_xz_yz_y = buffer_dddp[265];

    auto g_xz_xz_yz_z = buffer_dddp[266];

    auto g_xz_xz_zz_x = buffer_dddp[267];

    auto g_xz_xz_zz_y = buffer_dddp[268];

    auto g_xz_xz_zz_z = buffer_dddp[269];

    auto g_xz_yy_xx_x = buffer_dddp[270];

    auto g_xz_yy_xx_y = buffer_dddp[271];

    auto g_xz_yy_xx_z = buffer_dddp[272];

    auto g_xz_yy_xy_x = buffer_dddp[273];

    auto g_xz_yy_xy_y = buffer_dddp[274];

    auto g_xz_yy_xy_z = buffer_dddp[275];

    auto g_xz_yy_xz_x = buffer_dddp[276];

    auto g_xz_yy_xz_y = buffer_dddp[277];

    auto g_xz_yy_xz_z = buffer_dddp[278];

    auto g_xz_yy_yy_x = buffer_dddp[279];

    auto g_xz_yy_yy_y = buffer_dddp[280];

    auto g_xz_yy_yy_z = buffer_dddp[281];

    auto g_xz_yy_yz_x = buffer_dddp[282];

    auto g_xz_yy_yz_y = buffer_dddp[283];

    auto g_xz_yy_yz_z = buffer_dddp[284];

    auto g_xz_yy_zz_x = buffer_dddp[285];

    auto g_xz_yy_zz_y = buffer_dddp[286];

    auto g_xz_yy_zz_z = buffer_dddp[287];

    auto g_xz_yz_xx_x = buffer_dddp[288];

    auto g_xz_yz_xx_y = buffer_dddp[289];

    auto g_xz_yz_xx_z = buffer_dddp[290];

    auto g_xz_yz_xy_x = buffer_dddp[291];

    auto g_xz_yz_xy_y = buffer_dddp[292];

    auto g_xz_yz_xy_z = buffer_dddp[293];

    auto g_xz_yz_xz_x = buffer_dddp[294];

    auto g_xz_yz_xz_y = buffer_dddp[295];

    auto g_xz_yz_xz_z = buffer_dddp[296];

    auto g_xz_yz_yy_x = buffer_dddp[297];

    auto g_xz_yz_yy_y = buffer_dddp[298];

    auto g_xz_yz_yy_z = buffer_dddp[299];

    auto g_xz_yz_yz_x = buffer_dddp[300];

    auto g_xz_yz_yz_y = buffer_dddp[301];

    auto g_xz_yz_yz_z = buffer_dddp[302];

    auto g_xz_yz_zz_x = buffer_dddp[303];

    auto g_xz_yz_zz_y = buffer_dddp[304];

    auto g_xz_yz_zz_z = buffer_dddp[305];

    auto g_xz_zz_xx_x = buffer_dddp[306];

    auto g_xz_zz_xx_y = buffer_dddp[307];

    auto g_xz_zz_xx_z = buffer_dddp[308];

    auto g_xz_zz_xy_x = buffer_dddp[309];

    auto g_xz_zz_xy_y = buffer_dddp[310];

    auto g_xz_zz_xy_z = buffer_dddp[311];

    auto g_xz_zz_xz_x = buffer_dddp[312];

    auto g_xz_zz_xz_y = buffer_dddp[313];

    auto g_xz_zz_xz_z = buffer_dddp[314];

    auto g_xz_zz_yy_x = buffer_dddp[315];

    auto g_xz_zz_yy_y = buffer_dddp[316];

    auto g_xz_zz_yy_z = buffer_dddp[317];

    auto g_xz_zz_yz_x = buffer_dddp[318];

    auto g_xz_zz_yz_y = buffer_dddp[319];

    auto g_xz_zz_yz_z = buffer_dddp[320];

    auto g_xz_zz_zz_x = buffer_dddp[321];

    auto g_xz_zz_zz_y = buffer_dddp[322];

    auto g_xz_zz_zz_z = buffer_dddp[323];

    auto g_yy_xx_xx_x = buffer_dddp[324];

    auto g_yy_xx_xx_y = buffer_dddp[325];

    auto g_yy_xx_xx_z = buffer_dddp[326];

    auto g_yy_xx_xy_x = buffer_dddp[327];

    auto g_yy_xx_xy_y = buffer_dddp[328];

    auto g_yy_xx_xy_z = buffer_dddp[329];

    auto g_yy_xx_xz_x = buffer_dddp[330];

    auto g_yy_xx_xz_y = buffer_dddp[331];

    auto g_yy_xx_xz_z = buffer_dddp[332];

    auto g_yy_xx_yy_x = buffer_dddp[333];

    auto g_yy_xx_yy_y = buffer_dddp[334];

    auto g_yy_xx_yy_z = buffer_dddp[335];

    auto g_yy_xx_yz_x = buffer_dddp[336];

    auto g_yy_xx_yz_y = buffer_dddp[337];

    auto g_yy_xx_yz_z = buffer_dddp[338];

    auto g_yy_xx_zz_x = buffer_dddp[339];

    auto g_yy_xx_zz_y = buffer_dddp[340];

    auto g_yy_xx_zz_z = buffer_dddp[341];

    auto g_yy_xy_xx_x = buffer_dddp[342];

    auto g_yy_xy_xx_y = buffer_dddp[343];

    auto g_yy_xy_xx_z = buffer_dddp[344];

    auto g_yy_xy_xy_x = buffer_dddp[345];

    auto g_yy_xy_xy_y = buffer_dddp[346];

    auto g_yy_xy_xy_z = buffer_dddp[347];

    auto g_yy_xy_xz_x = buffer_dddp[348];

    auto g_yy_xy_xz_y = buffer_dddp[349];

    auto g_yy_xy_xz_z = buffer_dddp[350];

    auto g_yy_xy_yy_x = buffer_dddp[351];

    auto g_yy_xy_yy_y = buffer_dddp[352];

    auto g_yy_xy_yy_z = buffer_dddp[353];

    auto g_yy_xy_yz_x = buffer_dddp[354];

    auto g_yy_xy_yz_y = buffer_dddp[355];

    auto g_yy_xy_yz_z = buffer_dddp[356];

    auto g_yy_xy_zz_x = buffer_dddp[357];

    auto g_yy_xy_zz_y = buffer_dddp[358];

    auto g_yy_xy_zz_z = buffer_dddp[359];

    auto g_yy_xz_xx_x = buffer_dddp[360];

    auto g_yy_xz_xx_y = buffer_dddp[361];

    auto g_yy_xz_xx_z = buffer_dddp[362];

    auto g_yy_xz_xy_x = buffer_dddp[363];

    auto g_yy_xz_xy_y = buffer_dddp[364];

    auto g_yy_xz_xy_z = buffer_dddp[365];

    auto g_yy_xz_xz_x = buffer_dddp[366];

    auto g_yy_xz_xz_y = buffer_dddp[367];

    auto g_yy_xz_xz_z = buffer_dddp[368];

    auto g_yy_xz_yy_x = buffer_dddp[369];

    auto g_yy_xz_yy_y = buffer_dddp[370];

    auto g_yy_xz_yy_z = buffer_dddp[371];

    auto g_yy_xz_yz_x = buffer_dddp[372];

    auto g_yy_xz_yz_y = buffer_dddp[373];

    auto g_yy_xz_yz_z = buffer_dddp[374];

    auto g_yy_xz_zz_x = buffer_dddp[375];

    auto g_yy_xz_zz_y = buffer_dddp[376];

    auto g_yy_xz_zz_z = buffer_dddp[377];

    auto g_yy_yy_xx_x = buffer_dddp[378];

    auto g_yy_yy_xx_y = buffer_dddp[379];

    auto g_yy_yy_xx_z = buffer_dddp[380];

    auto g_yy_yy_xy_x = buffer_dddp[381];

    auto g_yy_yy_xy_y = buffer_dddp[382];

    auto g_yy_yy_xy_z = buffer_dddp[383];

    auto g_yy_yy_xz_x = buffer_dddp[384];

    auto g_yy_yy_xz_y = buffer_dddp[385];

    auto g_yy_yy_xz_z = buffer_dddp[386];

    auto g_yy_yy_yy_x = buffer_dddp[387];

    auto g_yy_yy_yy_y = buffer_dddp[388];

    auto g_yy_yy_yy_z = buffer_dddp[389];

    auto g_yy_yy_yz_x = buffer_dddp[390];

    auto g_yy_yy_yz_y = buffer_dddp[391];

    auto g_yy_yy_yz_z = buffer_dddp[392];

    auto g_yy_yy_zz_x = buffer_dddp[393];

    auto g_yy_yy_zz_y = buffer_dddp[394];

    auto g_yy_yy_zz_z = buffer_dddp[395];

    auto g_yy_yz_xx_x = buffer_dddp[396];

    auto g_yy_yz_xx_y = buffer_dddp[397];

    auto g_yy_yz_xx_z = buffer_dddp[398];

    auto g_yy_yz_xy_x = buffer_dddp[399];

    auto g_yy_yz_xy_y = buffer_dddp[400];

    auto g_yy_yz_xy_z = buffer_dddp[401];

    auto g_yy_yz_xz_x = buffer_dddp[402];

    auto g_yy_yz_xz_y = buffer_dddp[403];

    auto g_yy_yz_xz_z = buffer_dddp[404];

    auto g_yy_yz_yy_x = buffer_dddp[405];

    auto g_yy_yz_yy_y = buffer_dddp[406];

    auto g_yy_yz_yy_z = buffer_dddp[407];

    auto g_yy_yz_yz_x = buffer_dddp[408];

    auto g_yy_yz_yz_y = buffer_dddp[409];

    auto g_yy_yz_yz_z = buffer_dddp[410];

    auto g_yy_yz_zz_x = buffer_dddp[411];

    auto g_yy_yz_zz_y = buffer_dddp[412];

    auto g_yy_yz_zz_z = buffer_dddp[413];

    auto g_yy_zz_xx_x = buffer_dddp[414];

    auto g_yy_zz_xx_y = buffer_dddp[415];

    auto g_yy_zz_xx_z = buffer_dddp[416];

    auto g_yy_zz_xy_x = buffer_dddp[417];

    auto g_yy_zz_xy_y = buffer_dddp[418];

    auto g_yy_zz_xy_z = buffer_dddp[419];

    auto g_yy_zz_xz_x = buffer_dddp[420];

    auto g_yy_zz_xz_y = buffer_dddp[421];

    auto g_yy_zz_xz_z = buffer_dddp[422];

    auto g_yy_zz_yy_x = buffer_dddp[423];

    auto g_yy_zz_yy_y = buffer_dddp[424];

    auto g_yy_zz_yy_z = buffer_dddp[425];

    auto g_yy_zz_yz_x = buffer_dddp[426];

    auto g_yy_zz_yz_y = buffer_dddp[427];

    auto g_yy_zz_yz_z = buffer_dddp[428];

    auto g_yy_zz_zz_x = buffer_dddp[429];

    auto g_yy_zz_zz_y = buffer_dddp[430];

    auto g_yy_zz_zz_z = buffer_dddp[431];

    auto g_yz_xx_xx_x = buffer_dddp[432];

    auto g_yz_xx_xx_y = buffer_dddp[433];

    auto g_yz_xx_xx_z = buffer_dddp[434];

    auto g_yz_xx_xy_x = buffer_dddp[435];

    auto g_yz_xx_xy_y = buffer_dddp[436];

    auto g_yz_xx_xy_z = buffer_dddp[437];

    auto g_yz_xx_xz_x = buffer_dddp[438];

    auto g_yz_xx_xz_y = buffer_dddp[439];

    auto g_yz_xx_xz_z = buffer_dddp[440];

    auto g_yz_xx_yy_x = buffer_dddp[441];

    auto g_yz_xx_yy_y = buffer_dddp[442];

    auto g_yz_xx_yy_z = buffer_dddp[443];

    auto g_yz_xx_yz_x = buffer_dddp[444];

    auto g_yz_xx_yz_y = buffer_dddp[445];

    auto g_yz_xx_yz_z = buffer_dddp[446];

    auto g_yz_xx_zz_x = buffer_dddp[447];

    auto g_yz_xx_zz_y = buffer_dddp[448];

    auto g_yz_xx_zz_z = buffer_dddp[449];

    auto g_yz_xy_xx_x = buffer_dddp[450];

    auto g_yz_xy_xx_y = buffer_dddp[451];

    auto g_yz_xy_xx_z = buffer_dddp[452];

    auto g_yz_xy_xy_x = buffer_dddp[453];

    auto g_yz_xy_xy_y = buffer_dddp[454];

    auto g_yz_xy_xy_z = buffer_dddp[455];

    auto g_yz_xy_xz_x = buffer_dddp[456];

    auto g_yz_xy_xz_y = buffer_dddp[457];

    auto g_yz_xy_xz_z = buffer_dddp[458];

    auto g_yz_xy_yy_x = buffer_dddp[459];

    auto g_yz_xy_yy_y = buffer_dddp[460];

    auto g_yz_xy_yy_z = buffer_dddp[461];

    auto g_yz_xy_yz_x = buffer_dddp[462];

    auto g_yz_xy_yz_y = buffer_dddp[463];

    auto g_yz_xy_yz_z = buffer_dddp[464];

    auto g_yz_xy_zz_x = buffer_dddp[465];

    auto g_yz_xy_zz_y = buffer_dddp[466];

    auto g_yz_xy_zz_z = buffer_dddp[467];

    auto g_yz_xz_xx_x = buffer_dddp[468];

    auto g_yz_xz_xx_y = buffer_dddp[469];

    auto g_yz_xz_xx_z = buffer_dddp[470];

    auto g_yz_xz_xy_x = buffer_dddp[471];

    auto g_yz_xz_xy_y = buffer_dddp[472];

    auto g_yz_xz_xy_z = buffer_dddp[473];

    auto g_yz_xz_xz_x = buffer_dddp[474];

    auto g_yz_xz_xz_y = buffer_dddp[475];

    auto g_yz_xz_xz_z = buffer_dddp[476];

    auto g_yz_xz_yy_x = buffer_dddp[477];

    auto g_yz_xz_yy_y = buffer_dddp[478];

    auto g_yz_xz_yy_z = buffer_dddp[479];

    auto g_yz_xz_yz_x = buffer_dddp[480];

    auto g_yz_xz_yz_y = buffer_dddp[481];

    auto g_yz_xz_yz_z = buffer_dddp[482];

    auto g_yz_xz_zz_x = buffer_dddp[483];

    auto g_yz_xz_zz_y = buffer_dddp[484];

    auto g_yz_xz_zz_z = buffer_dddp[485];

    auto g_yz_yy_xx_x = buffer_dddp[486];

    auto g_yz_yy_xx_y = buffer_dddp[487];

    auto g_yz_yy_xx_z = buffer_dddp[488];

    auto g_yz_yy_xy_x = buffer_dddp[489];

    auto g_yz_yy_xy_y = buffer_dddp[490];

    auto g_yz_yy_xy_z = buffer_dddp[491];

    auto g_yz_yy_xz_x = buffer_dddp[492];

    auto g_yz_yy_xz_y = buffer_dddp[493];

    auto g_yz_yy_xz_z = buffer_dddp[494];

    auto g_yz_yy_yy_x = buffer_dddp[495];

    auto g_yz_yy_yy_y = buffer_dddp[496];

    auto g_yz_yy_yy_z = buffer_dddp[497];

    auto g_yz_yy_yz_x = buffer_dddp[498];

    auto g_yz_yy_yz_y = buffer_dddp[499];

    auto g_yz_yy_yz_z = buffer_dddp[500];

    auto g_yz_yy_zz_x = buffer_dddp[501];

    auto g_yz_yy_zz_y = buffer_dddp[502];

    auto g_yz_yy_zz_z = buffer_dddp[503];

    auto g_yz_yz_xx_x = buffer_dddp[504];

    auto g_yz_yz_xx_y = buffer_dddp[505];

    auto g_yz_yz_xx_z = buffer_dddp[506];

    auto g_yz_yz_xy_x = buffer_dddp[507];

    auto g_yz_yz_xy_y = buffer_dddp[508];

    auto g_yz_yz_xy_z = buffer_dddp[509];

    auto g_yz_yz_xz_x = buffer_dddp[510];

    auto g_yz_yz_xz_y = buffer_dddp[511];

    auto g_yz_yz_xz_z = buffer_dddp[512];

    auto g_yz_yz_yy_x = buffer_dddp[513];

    auto g_yz_yz_yy_y = buffer_dddp[514];

    auto g_yz_yz_yy_z = buffer_dddp[515];

    auto g_yz_yz_yz_x = buffer_dddp[516];

    auto g_yz_yz_yz_y = buffer_dddp[517];

    auto g_yz_yz_yz_z = buffer_dddp[518];

    auto g_yz_yz_zz_x = buffer_dddp[519];

    auto g_yz_yz_zz_y = buffer_dddp[520];

    auto g_yz_yz_zz_z = buffer_dddp[521];

    auto g_yz_zz_xx_x = buffer_dddp[522];

    auto g_yz_zz_xx_y = buffer_dddp[523];

    auto g_yz_zz_xx_z = buffer_dddp[524];

    auto g_yz_zz_xy_x = buffer_dddp[525];

    auto g_yz_zz_xy_y = buffer_dddp[526];

    auto g_yz_zz_xy_z = buffer_dddp[527];

    auto g_yz_zz_xz_x = buffer_dddp[528];

    auto g_yz_zz_xz_y = buffer_dddp[529];

    auto g_yz_zz_xz_z = buffer_dddp[530];

    auto g_yz_zz_yy_x = buffer_dddp[531];

    auto g_yz_zz_yy_y = buffer_dddp[532];

    auto g_yz_zz_yy_z = buffer_dddp[533];

    auto g_yz_zz_yz_x = buffer_dddp[534];

    auto g_yz_zz_yz_y = buffer_dddp[535];

    auto g_yz_zz_yz_z = buffer_dddp[536];

    auto g_yz_zz_zz_x = buffer_dddp[537];

    auto g_yz_zz_zz_y = buffer_dddp[538];

    auto g_yz_zz_zz_z = buffer_dddp[539];

    auto g_zz_xx_xx_x = buffer_dddp[540];

    auto g_zz_xx_xx_y = buffer_dddp[541];

    auto g_zz_xx_xx_z = buffer_dddp[542];

    auto g_zz_xx_xy_x = buffer_dddp[543];

    auto g_zz_xx_xy_y = buffer_dddp[544];

    auto g_zz_xx_xy_z = buffer_dddp[545];

    auto g_zz_xx_xz_x = buffer_dddp[546];

    auto g_zz_xx_xz_y = buffer_dddp[547];

    auto g_zz_xx_xz_z = buffer_dddp[548];

    auto g_zz_xx_yy_x = buffer_dddp[549];

    auto g_zz_xx_yy_y = buffer_dddp[550];

    auto g_zz_xx_yy_z = buffer_dddp[551];

    auto g_zz_xx_yz_x = buffer_dddp[552];

    auto g_zz_xx_yz_y = buffer_dddp[553];

    auto g_zz_xx_yz_z = buffer_dddp[554];

    auto g_zz_xx_zz_x = buffer_dddp[555];

    auto g_zz_xx_zz_y = buffer_dddp[556];

    auto g_zz_xx_zz_z = buffer_dddp[557];

    auto g_zz_xy_xx_x = buffer_dddp[558];

    auto g_zz_xy_xx_y = buffer_dddp[559];

    auto g_zz_xy_xx_z = buffer_dddp[560];

    auto g_zz_xy_xy_x = buffer_dddp[561];

    auto g_zz_xy_xy_y = buffer_dddp[562];

    auto g_zz_xy_xy_z = buffer_dddp[563];

    auto g_zz_xy_xz_x = buffer_dddp[564];

    auto g_zz_xy_xz_y = buffer_dddp[565];

    auto g_zz_xy_xz_z = buffer_dddp[566];

    auto g_zz_xy_yy_x = buffer_dddp[567];

    auto g_zz_xy_yy_y = buffer_dddp[568];

    auto g_zz_xy_yy_z = buffer_dddp[569];

    auto g_zz_xy_yz_x = buffer_dddp[570];

    auto g_zz_xy_yz_y = buffer_dddp[571];

    auto g_zz_xy_yz_z = buffer_dddp[572];

    auto g_zz_xy_zz_x = buffer_dddp[573];

    auto g_zz_xy_zz_y = buffer_dddp[574];

    auto g_zz_xy_zz_z = buffer_dddp[575];

    auto g_zz_xz_xx_x = buffer_dddp[576];

    auto g_zz_xz_xx_y = buffer_dddp[577];

    auto g_zz_xz_xx_z = buffer_dddp[578];

    auto g_zz_xz_xy_x = buffer_dddp[579];

    auto g_zz_xz_xy_y = buffer_dddp[580];

    auto g_zz_xz_xy_z = buffer_dddp[581];

    auto g_zz_xz_xz_x = buffer_dddp[582];

    auto g_zz_xz_xz_y = buffer_dddp[583];

    auto g_zz_xz_xz_z = buffer_dddp[584];

    auto g_zz_xz_yy_x = buffer_dddp[585];

    auto g_zz_xz_yy_y = buffer_dddp[586];

    auto g_zz_xz_yy_z = buffer_dddp[587];

    auto g_zz_xz_yz_x = buffer_dddp[588];

    auto g_zz_xz_yz_y = buffer_dddp[589];

    auto g_zz_xz_yz_z = buffer_dddp[590];

    auto g_zz_xz_zz_x = buffer_dddp[591];

    auto g_zz_xz_zz_y = buffer_dddp[592];

    auto g_zz_xz_zz_z = buffer_dddp[593];

    auto g_zz_yy_xx_x = buffer_dddp[594];

    auto g_zz_yy_xx_y = buffer_dddp[595];

    auto g_zz_yy_xx_z = buffer_dddp[596];

    auto g_zz_yy_xy_x = buffer_dddp[597];

    auto g_zz_yy_xy_y = buffer_dddp[598];

    auto g_zz_yy_xy_z = buffer_dddp[599];

    auto g_zz_yy_xz_x = buffer_dddp[600];

    auto g_zz_yy_xz_y = buffer_dddp[601];

    auto g_zz_yy_xz_z = buffer_dddp[602];

    auto g_zz_yy_yy_x = buffer_dddp[603];

    auto g_zz_yy_yy_y = buffer_dddp[604];

    auto g_zz_yy_yy_z = buffer_dddp[605];

    auto g_zz_yy_yz_x = buffer_dddp[606];

    auto g_zz_yy_yz_y = buffer_dddp[607];

    auto g_zz_yy_yz_z = buffer_dddp[608];

    auto g_zz_yy_zz_x = buffer_dddp[609];

    auto g_zz_yy_zz_y = buffer_dddp[610];

    auto g_zz_yy_zz_z = buffer_dddp[611];

    auto g_zz_yz_xx_x = buffer_dddp[612];

    auto g_zz_yz_xx_y = buffer_dddp[613];

    auto g_zz_yz_xx_z = buffer_dddp[614];

    auto g_zz_yz_xy_x = buffer_dddp[615];

    auto g_zz_yz_xy_y = buffer_dddp[616];

    auto g_zz_yz_xy_z = buffer_dddp[617];

    auto g_zz_yz_xz_x = buffer_dddp[618];

    auto g_zz_yz_xz_y = buffer_dddp[619];

    auto g_zz_yz_xz_z = buffer_dddp[620];

    auto g_zz_yz_yy_x = buffer_dddp[621];

    auto g_zz_yz_yy_y = buffer_dddp[622];

    auto g_zz_yz_yy_z = buffer_dddp[623];

    auto g_zz_yz_yz_x = buffer_dddp[624];

    auto g_zz_yz_yz_y = buffer_dddp[625];

    auto g_zz_yz_yz_z = buffer_dddp[626];

    auto g_zz_yz_zz_x = buffer_dddp[627];

    auto g_zz_yz_zz_y = buffer_dddp[628];

    auto g_zz_yz_zz_z = buffer_dddp[629];

    auto g_zz_zz_xx_x = buffer_dddp[630];

    auto g_zz_zz_xx_y = buffer_dddp[631];

    auto g_zz_zz_xx_z = buffer_dddp[632];

    auto g_zz_zz_xy_x = buffer_dddp[633];

    auto g_zz_zz_xy_y = buffer_dddp[634];

    auto g_zz_zz_xy_z = buffer_dddp[635];

    auto g_zz_zz_xz_x = buffer_dddp[636];

    auto g_zz_zz_xz_y = buffer_dddp[637];

    auto g_zz_zz_xz_z = buffer_dddp[638];

    auto g_zz_zz_yy_x = buffer_dddp[639];

    auto g_zz_zz_yy_y = buffer_dddp[640];

    auto g_zz_zz_yy_z = buffer_dddp[641];

    auto g_zz_zz_yz_x = buffer_dddp[642];

    auto g_zz_zz_yz_y = buffer_dddp[643];

    auto g_zz_zz_yz_z = buffer_dddp[644];

    auto g_zz_zz_zz_x = buffer_dddp[645];

    auto g_zz_zz_zz_y = buffer_dddp[646];

    auto g_zz_zz_zz_z = buffer_dddp[647];

    /// Set up components of integrals buffer : buffer_1010_pdpp

    auto g_x_0_x_0_x_xx_x_x = buffer_1010_pdpp[0];

    auto g_x_0_x_0_x_xx_x_y = buffer_1010_pdpp[1];

    auto g_x_0_x_0_x_xx_x_z = buffer_1010_pdpp[2];

    auto g_x_0_x_0_x_xx_y_x = buffer_1010_pdpp[3];

    auto g_x_0_x_0_x_xx_y_y = buffer_1010_pdpp[4];

    auto g_x_0_x_0_x_xx_y_z = buffer_1010_pdpp[5];

    auto g_x_0_x_0_x_xx_z_x = buffer_1010_pdpp[6];

    auto g_x_0_x_0_x_xx_z_y = buffer_1010_pdpp[7];

    auto g_x_0_x_0_x_xx_z_z = buffer_1010_pdpp[8];

    auto g_x_0_x_0_x_xy_x_x = buffer_1010_pdpp[9];

    auto g_x_0_x_0_x_xy_x_y = buffer_1010_pdpp[10];

    auto g_x_0_x_0_x_xy_x_z = buffer_1010_pdpp[11];

    auto g_x_0_x_0_x_xy_y_x = buffer_1010_pdpp[12];

    auto g_x_0_x_0_x_xy_y_y = buffer_1010_pdpp[13];

    auto g_x_0_x_0_x_xy_y_z = buffer_1010_pdpp[14];

    auto g_x_0_x_0_x_xy_z_x = buffer_1010_pdpp[15];

    auto g_x_0_x_0_x_xy_z_y = buffer_1010_pdpp[16];

    auto g_x_0_x_0_x_xy_z_z = buffer_1010_pdpp[17];

    auto g_x_0_x_0_x_xz_x_x = buffer_1010_pdpp[18];

    auto g_x_0_x_0_x_xz_x_y = buffer_1010_pdpp[19];

    auto g_x_0_x_0_x_xz_x_z = buffer_1010_pdpp[20];

    auto g_x_0_x_0_x_xz_y_x = buffer_1010_pdpp[21];

    auto g_x_0_x_0_x_xz_y_y = buffer_1010_pdpp[22];

    auto g_x_0_x_0_x_xz_y_z = buffer_1010_pdpp[23];

    auto g_x_0_x_0_x_xz_z_x = buffer_1010_pdpp[24];

    auto g_x_0_x_0_x_xz_z_y = buffer_1010_pdpp[25];

    auto g_x_0_x_0_x_xz_z_z = buffer_1010_pdpp[26];

    auto g_x_0_x_0_x_yy_x_x = buffer_1010_pdpp[27];

    auto g_x_0_x_0_x_yy_x_y = buffer_1010_pdpp[28];

    auto g_x_0_x_0_x_yy_x_z = buffer_1010_pdpp[29];

    auto g_x_0_x_0_x_yy_y_x = buffer_1010_pdpp[30];

    auto g_x_0_x_0_x_yy_y_y = buffer_1010_pdpp[31];

    auto g_x_0_x_0_x_yy_y_z = buffer_1010_pdpp[32];

    auto g_x_0_x_0_x_yy_z_x = buffer_1010_pdpp[33];

    auto g_x_0_x_0_x_yy_z_y = buffer_1010_pdpp[34];

    auto g_x_0_x_0_x_yy_z_z = buffer_1010_pdpp[35];

    auto g_x_0_x_0_x_yz_x_x = buffer_1010_pdpp[36];

    auto g_x_0_x_0_x_yz_x_y = buffer_1010_pdpp[37];

    auto g_x_0_x_0_x_yz_x_z = buffer_1010_pdpp[38];

    auto g_x_0_x_0_x_yz_y_x = buffer_1010_pdpp[39];

    auto g_x_0_x_0_x_yz_y_y = buffer_1010_pdpp[40];

    auto g_x_0_x_0_x_yz_y_z = buffer_1010_pdpp[41];

    auto g_x_0_x_0_x_yz_z_x = buffer_1010_pdpp[42];

    auto g_x_0_x_0_x_yz_z_y = buffer_1010_pdpp[43];

    auto g_x_0_x_0_x_yz_z_z = buffer_1010_pdpp[44];

    auto g_x_0_x_0_x_zz_x_x = buffer_1010_pdpp[45];

    auto g_x_0_x_0_x_zz_x_y = buffer_1010_pdpp[46];

    auto g_x_0_x_0_x_zz_x_z = buffer_1010_pdpp[47];

    auto g_x_0_x_0_x_zz_y_x = buffer_1010_pdpp[48];

    auto g_x_0_x_0_x_zz_y_y = buffer_1010_pdpp[49];

    auto g_x_0_x_0_x_zz_y_z = buffer_1010_pdpp[50];

    auto g_x_0_x_0_x_zz_z_x = buffer_1010_pdpp[51];

    auto g_x_0_x_0_x_zz_z_y = buffer_1010_pdpp[52];

    auto g_x_0_x_0_x_zz_z_z = buffer_1010_pdpp[53];

    auto g_x_0_x_0_y_xx_x_x = buffer_1010_pdpp[54];

    auto g_x_0_x_0_y_xx_x_y = buffer_1010_pdpp[55];

    auto g_x_0_x_0_y_xx_x_z = buffer_1010_pdpp[56];

    auto g_x_0_x_0_y_xx_y_x = buffer_1010_pdpp[57];

    auto g_x_0_x_0_y_xx_y_y = buffer_1010_pdpp[58];

    auto g_x_0_x_0_y_xx_y_z = buffer_1010_pdpp[59];

    auto g_x_0_x_0_y_xx_z_x = buffer_1010_pdpp[60];

    auto g_x_0_x_0_y_xx_z_y = buffer_1010_pdpp[61];

    auto g_x_0_x_0_y_xx_z_z = buffer_1010_pdpp[62];

    auto g_x_0_x_0_y_xy_x_x = buffer_1010_pdpp[63];

    auto g_x_0_x_0_y_xy_x_y = buffer_1010_pdpp[64];

    auto g_x_0_x_0_y_xy_x_z = buffer_1010_pdpp[65];

    auto g_x_0_x_0_y_xy_y_x = buffer_1010_pdpp[66];

    auto g_x_0_x_0_y_xy_y_y = buffer_1010_pdpp[67];

    auto g_x_0_x_0_y_xy_y_z = buffer_1010_pdpp[68];

    auto g_x_0_x_0_y_xy_z_x = buffer_1010_pdpp[69];

    auto g_x_0_x_0_y_xy_z_y = buffer_1010_pdpp[70];

    auto g_x_0_x_0_y_xy_z_z = buffer_1010_pdpp[71];

    auto g_x_0_x_0_y_xz_x_x = buffer_1010_pdpp[72];

    auto g_x_0_x_0_y_xz_x_y = buffer_1010_pdpp[73];

    auto g_x_0_x_0_y_xz_x_z = buffer_1010_pdpp[74];

    auto g_x_0_x_0_y_xz_y_x = buffer_1010_pdpp[75];

    auto g_x_0_x_0_y_xz_y_y = buffer_1010_pdpp[76];

    auto g_x_0_x_0_y_xz_y_z = buffer_1010_pdpp[77];

    auto g_x_0_x_0_y_xz_z_x = buffer_1010_pdpp[78];

    auto g_x_0_x_0_y_xz_z_y = buffer_1010_pdpp[79];

    auto g_x_0_x_0_y_xz_z_z = buffer_1010_pdpp[80];

    auto g_x_0_x_0_y_yy_x_x = buffer_1010_pdpp[81];

    auto g_x_0_x_0_y_yy_x_y = buffer_1010_pdpp[82];

    auto g_x_0_x_0_y_yy_x_z = buffer_1010_pdpp[83];

    auto g_x_0_x_0_y_yy_y_x = buffer_1010_pdpp[84];

    auto g_x_0_x_0_y_yy_y_y = buffer_1010_pdpp[85];

    auto g_x_0_x_0_y_yy_y_z = buffer_1010_pdpp[86];

    auto g_x_0_x_0_y_yy_z_x = buffer_1010_pdpp[87];

    auto g_x_0_x_0_y_yy_z_y = buffer_1010_pdpp[88];

    auto g_x_0_x_0_y_yy_z_z = buffer_1010_pdpp[89];

    auto g_x_0_x_0_y_yz_x_x = buffer_1010_pdpp[90];

    auto g_x_0_x_0_y_yz_x_y = buffer_1010_pdpp[91];

    auto g_x_0_x_0_y_yz_x_z = buffer_1010_pdpp[92];

    auto g_x_0_x_0_y_yz_y_x = buffer_1010_pdpp[93];

    auto g_x_0_x_0_y_yz_y_y = buffer_1010_pdpp[94];

    auto g_x_0_x_0_y_yz_y_z = buffer_1010_pdpp[95];

    auto g_x_0_x_0_y_yz_z_x = buffer_1010_pdpp[96];

    auto g_x_0_x_0_y_yz_z_y = buffer_1010_pdpp[97];

    auto g_x_0_x_0_y_yz_z_z = buffer_1010_pdpp[98];

    auto g_x_0_x_0_y_zz_x_x = buffer_1010_pdpp[99];

    auto g_x_0_x_0_y_zz_x_y = buffer_1010_pdpp[100];

    auto g_x_0_x_0_y_zz_x_z = buffer_1010_pdpp[101];

    auto g_x_0_x_0_y_zz_y_x = buffer_1010_pdpp[102];

    auto g_x_0_x_0_y_zz_y_y = buffer_1010_pdpp[103];

    auto g_x_0_x_0_y_zz_y_z = buffer_1010_pdpp[104];

    auto g_x_0_x_0_y_zz_z_x = buffer_1010_pdpp[105];

    auto g_x_0_x_0_y_zz_z_y = buffer_1010_pdpp[106];

    auto g_x_0_x_0_y_zz_z_z = buffer_1010_pdpp[107];

    auto g_x_0_x_0_z_xx_x_x = buffer_1010_pdpp[108];

    auto g_x_0_x_0_z_xx_x_y = buffer_1010_pdpp[109];

    auto g_x_0_x_0_z_xx_x_z = buffer_1010_pdpp[110];

    auto g_x_0_x_0_z_xx_y_x = buffer_1010_pdpp[111];

    auto g_x_0_x_0_z_xx_y_y = buffer_1010_pdpp[112];

    auto g_x_0_x_0_z_xx_y_z = buffer_1010_pdpp[113];

    auto g_x_0_x_0_z_xx_z_x = buffer_1010_pdpp[114];

    auto g_x_0_x_0_z_xx_z_y = buffer_1010_pdpp[115];

    auto g_x_0_x_0_z_xx_z_z = buffer_1010_pdpp[116];

    auto g_x_0_x_0_z_xy_x_x = buffer_1010_pdpp[117];

    auto g_x_0_x_0_z_xy_x_y = buffer_1010_pdpp[118];

    auto g_x_0_x_0_z_xy_x_z = buffer_1010_pdpp[119];

    auto g_x_0_x_0_z_xy_y_x = buffer_1010_pdpp[120];

    auto g_x_0_x_0_z_xy_y_y = buffer_1010_pdpp[121];

    auto g_x_0_x_0_z_xy_y_z = buffer_1010_pdpp[122];

    auto g_x_0_x_0_z_xy_z_x = buffer_1010_pdpp[123];

    auto g_x_0_x_0_z_xy_z_y = buffer_1010_pdpp[124];

    auto g_x_0_x_0_z_xy_z_z = buffer_1010_pdpp[125];

    auto g_x_0_x_0_z_xz_x_x = buffer_1010_pdpp[126];

    auto g_x_0_x_0_z_xz_x_y = buffer_1010_pdpp[127];

    auto g_x_0_x_0_z_xz_x_z = buffer_1010_pdpp[128];

    auto g_x_0_x_0_z_xz_y_x = buffer_1010_pdpp[129];

    auto g_x_0_x_0_z_xz_y_y = buffer_1010_pdpp[130];

    auto g_x_0_x_0_z_xz_y_z = buffer_1010_pdpp[131];

    auto g_x_0_x_0_z_xz_z_x = buffer_1010_pdpp[132];

    auto g_x_0_x_0_z_xz_z_y = buffer_1010_pdpp[133];

    auto g_x_0_x_0_z_xz_z_z = buffer_1010_pdpp[134];

    auto g_x_0_x_0_z_yy_x_x = buffer_1010_pdpp[135];

    auto g_x_0_x_0_z_yy_x_y = buffer_1010_pdpp[136];

    auto g_x_0_x_0_z_yy_x_z = buffer_1010_pdpp[137];

    auto g_x_0_x_0_z_yy_y_x = buffer_1010_pdpp[138];

    auto g_x_0_x_0_z_yy_y_y = buffer_1010_pdpp[139];

    auto g_x_0_x_0_z_yy_y_z = buffer_1010_pdpp[140];

    auto g_x_0_x_0_z_yy_z_x = buffer_1010_pdpp[141];

    auto g_x_0_x_0_z_yy_z_y = buffer_1010_pdpp[142];

    auto g_x_0_x_0_z_yy_z_z = buffer_1010_pdpp[143];

    auto g_x_0_x_0_z_yz_x_x = buffer_1010_pdpp[144];

    auto g_x_0_x_0_z_yz_x_y = buffer_1010_pdpp[145];

    auto g_x_0_x_0_z_yz_x_z = buffer_1010_pdpp[146];

    auto g_x_0_x_0_z_yz_y_x = buffer_1010_pdpp[147];

    auto g_x_0_x_0_z_yz_y_y = buffer_1010_pdpp[148];

    auto g_x_0_x_0_z_yz_y_z = buffer_1010_pdpp[149];

    auto g_x_0_x_0_z_yz_z_x = buffer_1010_pdpp[150];

    auto g_x_0_x_0_z_yz_z_y = buffer_1010_pdpp[151];

    auto g_x_0_x_0_z_yz_z_z = buffer_1010_pdpp[152];

    auto g_x_0_x_0_z_zz_x_x = buffer_1010_pdpp[153];

    auto g_x_0_x_0_z_zz_x_y = buffer_1010_pdpp[154];

    auto g_x_0_x_0_z_zz_x_z = buffer_1010_pdpp[155];

    auto g_x_0_x_0_z_zz_y_x = buffer_1010_pdpp[156];

    auto g_x_0_x_0_z_zz_y_y = buffer_1010_pdpp[157];

    auto g_x_0_x_0_z_zz_y_z = buffer_1010_pdpp[158];

    auto g_x_0_x_0_z_zz_z_x = buffer_1010_pdpp[159];

    auto g_x_0_x_0_z_zz_z_y = buffer_1010_pdpp[160];

    auto g_x_0_x_0_z_zz_z_z = buffer_1010_pdpp[161];

    auto g_x_0_y_0_x_xx_x_x = buffer_1010_pdpp[162];

    auto g_x_0_y_0_x_xx_x_y = buffer_1010_pdpp[163];

    auto g_x_0_y_0_x_xx_x_z = buffer_1010_pdpp[164];

    auto g_x_0_y_0_x_xx_y_x = buffer_1010_pdpp[165];

    auto g_x_0_y_0_x_xx_y_y = buffer_1010_pdpp[166];

    auto g_x_0_y_0_x_xx_y_z = buffer_1010_pdpp[167];

    auto g_x_0_y_0_x_xx_z_x = buffer_1010_pdpp[168];

    auto g_x_0_y_0_x_xx_z_y = buffer_1010_pdpp[169];

    auto g_x_0_y_0_x_xx_z_z = buffer_1010_pdpp[170];

    auto g_x_0_y_0_x_xy_x_x = buffer_1010_pdpp[171];

    auto g_x_0_y_0_x_xy_x_y = buffer_1010_pdpp[172];

    auto g_x_0_y_0_x_xy_x_z = buffer_1010_pdpp[173];

    auto g_x_0_y_0_x_xy_y_x = buffer_1010_pdpp[174];

    auto g_x_0_y_0_x_xy_y_y = buffer_1010_pdpp[175];

    auto g_x_0_y_0_x_xy_y_z = buffer_1010_pdpp[176];

    auto g_x_0_y_0_x_xy_z_x = buffer_1010_pdpp[177];

    auto g_x_0_y_0_x_xy_z_y = buffer_1010_pdpp[178];

    auto g_x_0_y_0_x_xy_z_z = buffer_1010_pdpp[179];

    auto g_x_0_y_0_x_xz_x_x = buffer_1010_pdpp[180];

    auto g_x_0_y_0_x_xz_x_y = buffer_1010_pdpp[181];

    auto g_x_0_y_0_x_xz_x_z = buffer_1010_pdpp[182];

    auto g_x_0_y_0_x_xz_y_x = buffer_1010_pdpp[183];

    auto g_x_0_y_0_x_xz_y_y = buffer_1010_pdpp[184];

    auto g_x_0_y_0_x_xz_y_z = buffer_1010_pdpp[185];

    auto g_x_0_y_0_x_xz_z_x = buffer_1010_pdpp[186];

    auto g_x_0_y_0_x_xz_z_y = buffer_1010_pdpp[187];

    auto g_x_0_y_0_x_xz_z_z = buffer_1010_pdpp[188];

    auto g_x_0_y_0_x_yy_x_x = buffer_1010_pdpp[189];

    auto g_x_0_y_0_x_yy_x_y = buffer_1010_pdpp[190];

    auto g_x_0_y_0_x_yy_x_z = buffer_1010_pdpp[191];

    auto g_x_0_y_0_x_yy_y_x = buffer_1010_pdpp[192];

    auto g_x_0_y_0_x_yy_y_y = buffer_1010_pdpp[193];

    auto g_x_0_y_0_x_yy_y_z = buffer_1010_pdpp[194];

    auto g_x_0_y_0_x_yy_z_x = buffer_1010_pdpp[195];

    auto g_x_0_y_0_x_yy_z_y = buffer_1010_pdpp[196];

    auto g_x_0_y_0_x_yy_z_z = buffer_1010_pdpp[197];

    auto g_x_0_y_0_x_yz_x_x = buffer_1010_pdpp[198];

    auto g_x_0_y_0_x_yz_x_y = buffer_1010_pdpp[199];

    auto g_x_0_y_0_x_yz_x_z = buffer_1010_pdpp[200];

    auto g_x_0_y_0_x_yz_y_x = buffer_1010_pdpp[201];

    auto g_x_0_y_0_x_yz_y_y = buffer_1010_pdpp[202];

    auto g_x_0_y_0_x_yz_y_z = buffer_1010_pdpp[203];

    auto g_x_0_y_0_x_yz_z_x = buffer_1010_pdpp[204];

    auto g_x_0_y_0_x_yz_z_y = buffer_1010_pdpp[205];

    auto g_x_0_y_0_x_yz_z_z = buffer_1010_pdpp[206];

    auto g_x_0_y_0_x_zz_x_x = buffer_1010_pdpp[207];

    auto g_x_0_y_0_x_zz_x_y = buffer_1010_pdpp[208];

    auto g_x_0_y_0_x_zz_x_z = buffer_1010_pdpp[209];

    auto g_x_0_y_0_x_zz_y_x = buffer_1010_pdpp[210];

    auto g_x_0_y_0_x_zz_y_y = buffer_1010_pdpp[211];

    auto g_x_0_y_0_x_zz_y_z = buffer_1010_pdpp[212];

    auto g_x_0_y_0_x_zz_z_x = buffer_1010_pdpp[213];

    auto g_x_0_y_0_x_zz_z_y = buffer_1010_pdpp[214];

    auto g_x_0_y_0_x_zz_z_z = buffer_1010_pdpp[215];

    auto g_x_0_y_0_y_xx_x_x = buffer_1010_pdpp[216];

    auto g_x_0_y_0_y_xx_x_y = buffer_1010_pdpp[217];

    auto g_x_0_y_0_y_xx_x_z = buffer_1010_pdpp[218];

    auto g_x_0_y_0_y_xx_y_x = buffer_1010_pdpp[219];

    auto g_x_0_y_0_y_xx_y_y = buffer_1010_pdpp[220];

    auto g_x_0_y_0_y_xx_y_z = buffer_1010_pdpp[221];

    auto g_x_0_y_0_y_xx_z_x = buffer_1010_pdpp[222];

    auto g_x_0_y_0_y_xx_z_y = buffer_1010_pdpp[223];

    auto g_x_0_y_0_y_xx_z_z = buffer_1010_pdpp[224];

    auto g_x_0_y_0_y_xy_x_x = buffer_1010_pdpp[225];

    auto g_x_0_y_0_y_xy_x_y = buffer_1010_pdpp[226];

    auto g_x_0_y_0_y_xy_x_z = buffer_1010_pdpp[227];

    auto g_x_0_y_0_y_xy_y_x = buffer_1010_pdpp[228];

    auto g_x_0_y_0_y_xy_y_y = buffer_1010_pdpp[229];

    auto g_x_0_y_0_y_xy_y_z = buffer_1010_pdpp[230];

    auto g_x_0_y_0_y_xy_z_x = buffer_1010_pdpp[231];

    auto g_x_0_y_0_y_xy_z_y = buffer_1010_pdpp[232];

    auto g_x_0_y_0_y_xy_z_z = buffer_1010_pdpp[233];

    auto g_x_0_y_0_y_xz_x_x = buffer_1010_pdpp[234];

    auto g_x_0_y_0_y_xz_x_y = buffer_1010_pdpp[235];

    auto g_x_0_y_0_y_xz_x_z = buffer_1010_pdpp[236];

    auto g_x_0_y_0_y_xz_y_x = buffer_1010_pdpp[237];

    auto g_x_0_y_0_y_xz_y_y = buffer_1010_pdpp[238];

    auto g_x_0_y_0_y_xz_y_z = buffer_1010_pdpp[239];

    auto g_x_0_y_0_y_xz_z_x = buffer_1010_pdpp[240];

    auto g_x_0_y_0_y_xz_z_y = buffer_1010_pdpp[241];

    auto g_x_0_y_0_y_xz_z_z = buffer_1010_pdpp[242];

    auto g_x_0_y_0_y_yy_x_x = buffer_1010_pdpp[243];

    auto g_x_0_y_0_y_yy_x_y = buffer_1010_pdpp[244];

    auto g_x_0_y_0_y_yy_x_z = buffer_1010_pdpp[245];

    auto g_x_0_y_0_y_yy_y_x = buffer_1010_pdpp[246];

    auto g_x_0_y_0_y_yy_y_y = buffer_1010_pdpp[247];

    auto g_x_0_y_0_y_yy_y_z = buffer_1010_pdpp[248];

    auto g_x_0_y_0_y_yy_z_x = buffer_1010_pdpp[249];

    auto g_x_0_y_0_y_yy_z_y = buffer_1010_pdpp[250];

    auto g_x_0_y_0_y_yy_z_z = buffer_1010_pdpp[251];

    auto g_x_0_y_0_y_yz_x_x = buffer_1010_pdpp[252];

    auto g_x_0_y_0_y_yz_x_y = buffer_1010_pdpp[253];

    auto g_x_0_y_0_y_yz_x_z = buffer_1010_pdpp[254];

    auto g_x_0_y_0_y_yz_y_x = buffer_1010_pdpp[255];

    auto g_x_0_y_0_y_yz_y_y = buffer_1010_pdpp[256];

    auto g_x_0_y_0_y_yz_y_z = buffer_1010_pdpp[257];

    auto g_x_0_y_0_y_yz_z_x = buffer_1010_pdpp[258];

    auto g_x_0_y_0_y_yz_z_y = buffer_1010_pdpp[259];

    auto g_x_0_y_0_y_yz_z_z = buffer_1010_pdpp[260];

    auto g_x_0_y_0_y_zz_x_x = buffer_1010_pdpp[261];

    auto g_x_0_y_0_y_zz_x_y = buffer_1010_pdpp[262];

    auto g_x_0_y_0_y_zz_x_z = buffer_1010_pdpp[263];

    auto g_x_0_y_0_y_zz_y_x = buffer_1010_pdpp[264];

    auto g_x_0_y_0_y_zz_y_y = buffer_1010_pdpp[265];

    auto g_x_0_y_0_y_zz_y_z = buffer_1010_pdpp[266];

    auto g_x_0_y_0_y_zz_z_x = buffer_1010_pdpp[267];

    auto g_x_0_y_0_y_zz_z_y = buffer_1010_pdpp[268];

    auto g_x_0_y_0_y_zz_z_z = buffer_1010_pdpp[269];

    auto g_x_0_y_0_z_xx_x_x = buffer_1010_pdpp[270];

    auto g_x_0_y_0_z_xx_x_y = buffer_1010_pdpp[271];

    auto g_x_0_y_0_z_xx_x_z = buffer_1010_pdpp[272];

    auto g_x_0_y_0_z_xx_y_x = buffer_1010_pdpp[273];

    auto g_x_0_y_0_z_xx_y_y = buffer_1010_pdpp[274];

    auto g_x_0_y_0_z_xx_y_z = buffer_1010_pdpp[275];

    auto g_x_0_y_0_z_xx_z_x = buffer_1010_pdpp[276];

    auto g_x_0_y_0_z_xx_z_y = buffer_1010_pdpp[277];

    auto g_x_0_y_0_z_xx_z_z = buffer_1010_pdpp[278];

    auto g_x_0_y_0_z_xy_x_x = buffer_1010_pdpp[279];

    auto g_x_0_y_0_z_xy_x_y = buffer_1010_pdpp[280];

    auto g_x_0_y_0_z_xy_x_z = buffer_1010_pdpp[281];

    auto g_x_0_y_0_z_xy_y_x = buffer_1010_pdpp[282];

    auto g_x_0_y_0_z_xy_y_y = buffer_1010_pdpp[283];

    auto g_x_0_y_0_z_xy_y_z = buffer_1010_pdpp[284];

    auto g_x_0_y_0_z_xy_z_x = buffer_1010_pdpp[285];

    auto g_x_0_y_0_z_xy_z_y = buffer_1010_pdpp[286];

    auto g_x_0_y_0_z_xy_z_z = buffer_1010_pdpp[287];

    auto g_x_0_y_0_z_xz_x_x = buffer_1010_pdpp[288];

    auto g_x_0_y_0_z_xz_x_y = buffer_1010_pdpp[289];

    auto g_x_0_y_0_z_xz_x_z = buffer_1010_pdpp[290];

    auto g_x_0_y_0_z_xz_y_x = buffer_1010_pdpp[291];

    auto g_x_0_y_0_z_xz_y_y = buffer_1010_pdpp[292];

    auto g_x_0_y_0_z_xz_y_z = buffer_1010_pdpp[293];

    auto g_x_0_y_0_z_xz_z_x = buffer_1010_pdpp[294];

    auto g_x_0_y_0_z_xz_z_y = buffer_1010_pdpp[295];

    auto g_x_0_y_0_z_xz_z_z = buffer_1010_pdpp[296];

    auto g_x_0_y_0_z_yy_x_x = buffer_1010_pdpp[297];

    auto g_x_0_y_0_z_yy_x_y = buffer_1010_pdpp[298];

    auto g_x_0_y_0_z_yy_x_z = buffer_1010_pdpp[299];

    auto g_x_0_y_0_z_yy_y_x = buffer_1010_pdpp[300];

    auto g_x_0_y_0_z_yy_y_y = buffer_1010_pdpp[301];

    auto g_x_0_y_0_z_yy_y_z = buffer_1010_pdpp[302];

    auto g_x_0_y_0_z_yy_z_x = buffer_1010_pdpp[303];

    auto g_x_0_y_0_z_yy_z_y = buffer_1010_pdpp[304];

    auto g_x_0_y_0_z_yy_z_z = buffer_1010_pdpp[305];

    auto g_x_0_y_0_z_yz_x_x = buffer_1010_pdpp[306];

    auto g_x_0_y_0_z_yz_x_y = buffer_1010_pdpp[307];

    auto g_x_0_y_0_z_yz_x_z = buffer_1010_pdpp[308];

    auto g_x_0_y_0_z_yz_y_x = buffer_1010_pdpp[309];

    auto g_x_0_y_0_z_yz_y_y = buffer_1010_pdpp[310];

    auto g_x_0_y_0_z_yz_y_z = buffer_1010_pdpp[311];

    auto g_x_0_y_0_z_yz_z_x = buffer_1010_pdpp[312];

    auto g_x_0_y_0_z_yz_z_y = buffer_1010_pdpp[313];

    auto g_x_0_y_0_z_yz_z_z = buffer_1010_pdpp[314];

    auto g_x_0_y_0_z_zz_x_x = buffer_1010_pdpp[315];

    auto g_x_0_y_0_z_zz_x_y = buffer_1010_pdpp[316];

    auto g_x_0_y_0_z_zz_x_z = buffer_1010_pdpp[317];

    auto g_x_0_y_0_z_zz_y_x = buffer_1010_pdpp[318];

    auto g_x_0_y_0_z_zz_y_y = buffer_1010_pdpp[319];

    auto g_x_0_y_0_z_zz_y_z = buffer_1010_pdpp[320];

    auto g_x_0_y_0_z_zz_z_x = buffer_1010_pdpp[321];

    auto g_x_0_y_0_z_zz_z_y = buffer_1010_pdpp[322];

    auto g_x_0_y_0_z_zz_z_z = buffer_1010_pdpp[323];

    auto g_x_0_z_0_x_xx_x_x = buffer_1010_pdpp[324];

    auto g_x_0_z_0_x_xx_x_y = buffer_1010_pdpp[325];

    auto g_x_0_z_0_x_xx_x_z = buffer_1010_pdpp[326];

    auto g_x_0_z_0_x_xx_y_x = buffer_1010_pdpp[327];

    auto g_x_0_z_0_x_xx_y_y = buffer_1010_pdpp[328];

    auto g_x_0_z_0_x_xx_y_z = buffer_1010_pdpp[329];

    auto g_x_0_z_0_x_xx_z_x = buffer_1010_pdpp[330];

    auto g_x_0_z_0_x_xx_z_y = buffer_1010_pdpp[331];

    auto g_x_0_z_0_x_xx_z_z = buffer_1010_pdpp[332];

    auto g_x_0_z_0_x_xy_x_x = buffer_1010_pdpp[333];

    auto g_x_0_z_0_x_xy_x_y = buffer_1010_pdpp[334];

    auto g_x_0_z_0_x_xy_x_z = buffer_1010_pdpp[335];

    auto g_x_0_z_0_x_xy_y_x = buffer_1010_pdpp[336];

    auto g_x_0_z_0_x_xy_y_y = buffer_1010_pdpp[337];

    auto g_x_0_z_0_x_xy_y_z = buffer_1010_pdpp[338];

    auto g_x_0_z_0_x_xy_z_x = buffer_1010_pdpp[339];

    auto g_x_0_z_0_x_xy_z_y = buffer_1010_pdpp[340];

    auto g_x_0_z_0_x_xy_z_z = buffer_1010_pdpp[341];

    auto g_x_0_z_0_x_xz_x_x = buffer_1010_pdpp[342];

    auto g_x_0_z_0_x_xz_x_y = buffer_1010_pdpp[343];

    auto g_x_0_z_0_x_xz_x_z = buffer_1010_pdpp[344];

    auto g_x_0_z_0_x_xz_y_x = buffer_1010_pdpp[345];

    auto g_x_0_z_0_x_xz_y_y = buffer_1010_pdpp[346];

    auto g_x_0_z_0_x_xz_y_z = buffer_1010_pdpp[347];

    auto g_x_0_z_0_x_xz_z_x = buffer_1010_pdpp[348];

    auto g_x_0_z_0_x_xz_z_y = buffer_1010_pdpp[349];

    auto g_x_0_z_0_x_xz_z_z = buffer_1010_pdpp[350];

    auto g_x_0_z_0_x_yy_x_x = buffer_1010_pdpp[351];

    auto g_x_0_z_0_x_yy_x_y = buffer_1010_pdpp[352];

    auto g_x_0_z_0_x_yy_x_z = buffer_1010_pdpp[353];

    auto g_x_0_z_0_x_yy_y_x = buffer_1010_pdpp[354];

    auto g_x_0_z_0_x_yy_y_y = buffer_1010_pdpp[355];

    auto g_x_0_z_0_x_yy_y_z = buffer_1010_pdpp[356];

    auto g_x_0_z_0_x_yy_z_x = buffer_1010_pdpp[357];

    auto g_x_0_z_0_x_yy_z_y = buffer_1010_pdpp[358];

    auto g_x_0_z_0_x_yy_z_z = buffer_1010_pdpp[359];

    auto g_x_0_z_0_x_yz_x_x = buffer_1010_pdpp[360];

    auto g_x_0_z_0_x_yz_x_y = buffer_1010_pdpp[361];

    auto g_x_0_z_0_x_yz_x_z = buffer_1010_pdpp[362];

    auto g_x_0_z_0_x_yz_y_x = buffer_1010_pdpp[363];

    auto g_x_0_z_0_x_yz_y_y = buffer_1010_pdpp[364];

    auto g_x_0_z_0_x_yz_y_z = buffer_1010_pdpp[365];

    auto g_x_0_z_0_x_yz_z_x = buffer_1010_pdpp[366];

    auto g_x_0_z_0_x_yz_z_y = buffer_1010_pdpp[367];

    auto g_x_0_z_0_x_yz_z_z = buffer_1010_pdpp[368];

    auto g_x_0_z_0_x_zz_x_x = buffer_1010_pdpp[369];

    auto g_x_0_z_0_x_zz_x_y = buffer_1010_pdpp[370];

    auto g_x_0_z_0_x_zz_x_z = buffer_1010_pdpp[371];

    auto g_x_0_z_0_x_zz_y_x = buffer_1010_pdpp[372];

    auto g_x_0_z_0_x_zz_y_y = buffer_1010_pdpp[373];

    auto g_x_0_z_0_x_zz_y_z = buffer_1010_pdpp[374];

    auto g_x_0_z_0_x_zz_z_x = buffer_1010_pdpp[375];

    auto g_x_0_z_0_x_zz_z_y = buffer_1010_pdpp[376];

    auto g_x_0_z_0_x_zz_z_z = buffer_1010_pdpp[377];

    auto g_x_0_z_0_y_xx_x_x = buffer_1010_pdpp[378];

    auto g_x_0_z_0_y_xx_x_y = buffer_1010_pdpp[379];

    auto g_x_0_z_0_y_xx_x_z = buffer_1010_pdpp[380];

    auto g_x_0_z_0_y_xx_y_x = buffer_1010_pdpp[381];

    auto g_x_0_z_0_y_xx_y_y = buffer_1010_pdpp[382];

    auto g_x_0_z_0_y_xx_y_z = buffer_1010_pdpp[383];

    auto g_x_0_z_0_y_xx_z_x = buffer_1010_pdpp[384];

    auto g_x_0_z_0_y_xx_z_y = buffer_1010_pdpp[385];

    auto g_x_0_z_0_y_xx_z_z = buffer_1010_pdpp[386];

    auto g_x_0_z_0_y_xy_x_x = buffer_1010_pdpp[387];

    auto g_x_0_z_0_y_xy_x_y = buffer_1010_pdpp[388];

    auto g_x_0_z_0_y_xy_x_z = buffer_1010_pdpp[389];

    auto g_x_0_z_0_y_xy_y_x = buffer_1010_pdpp[390];

    auto g_x_0_z_0_y_xy_y_y = buffer_1010_pdpp[391];

    auto g_x_0_z_0_y_xy_y_z = buffer_1010_pdpp[392];

    auto g_x_0_z_0_y_xy_z_x = buffer_1010_pdpp[393];

    auto g_x_0_z_0_y_xy_z_y = buffer_1010_pdpp[394];

    auto g_x_0_z_0_y_xy_z_z = buffer_1010_pdpp[395];

    auto g_x_0_z_0_y_xz_x_x = buffer_1010_pdpp[396];

    auto g_x_0_z_0_y_xz_x_y = buffer_1010_pdpp[397];

    auto g_x_0_z_0_y_xz_x_z = buffer_1010_pdpp[398];

    auto g_x_0_z_0_y_xz_y_x = buffer_1010_pdpp[399];

    auto g_x_0_z_0_y_xz_y_y = buffer_1010_pdpp[400];

    auto g_x_0_z_0_y_xz_y_z = buffer_1010_pdpp[401];

    auto g_x_0_z_0_y_xz_z_x = buffer_1010_pdpp[402];

    auto g_x_0_z_0_y_xz_z_y = buffer_1010_pdpp[403];

    auto g_x_0_z_0_y_xz_z_z = buffer_1010_pdpp[404];

    auto g_x_0_z_0_y_yy_x_x = buffer_1010_pdpp[405];

    auto g_x_0_z_0_y_yy_x_y = buffer_1010_pdpp[406];

    auto g_x_0_z_0_y_yy_x_z = buffer_1010_pdpp[407];

    auto g_x_0_z_0_y_yy_y_x = buffer_1010_pdpp[408];

    auto g_x_0_z_0_y_yy_y_y = buffer_1010_pdpp[409];

    auto g_x_0_z_0_y_yy_y_z = buffer_1010_pdpp[410];

    auto g_x_0_z_0_y_yy_z_x = buffer_1010_pdpp[411];

    auto g_x_0_z_0_y_yy_z_y = buffer_1010_pdpp[412];

    auto g_x_0_z_0_y_yy_z_z = buffer_1010_pdpp[413];

    auto g_x_0_z_0_y_yz_x_x = buffer_1010_pdpp[414];

    auto g_x_0_z_0_y_yz_x_y = buffer_1010_pdpp[415];

    auto g_x_0_z_0_y_yz_x_z = buffer_1010_pdpp[416];

    auto g_x_0_z_0_y_yz_y_x = buffer_1010_pdpp[417];

    auto g_x_0_z_0_y_yz_y_y = buffer_1010_pdpp[418];

    auto g_x_0_z_0_y_yz_y_z = buffer_1010_pdpp[419];

    auto g_x_0_z_0_y_yz_z_x = buffer_1010_pdpp[420];

    auto g_x_0_z_0_y_yz_z_y = buffer_1010_pdpp[421];

    auto g_x_0_z_0_y_yz_z_z = buffer_1010_pdpp[422];

    auto g_x_0_z_0_y_zz_x_x = buffer_1010_pdpp[423];

    auto g_x_0_z_0_y_zz_x_y = buffer_1010_pdpp[424];

    auto g_x_0_z_0_y_zz_x_z = buffer_1010_pdpp[425];

    auto g_x_0_z_0_y_zz_y_x = buffer_1010_pdpp[426];

    auto g_x_0_z_0_y_zz_y_y = buffer_1010_pdpp[427];

    auto g_x_0_z_0_y_zz_y_z = buffer_1010_pdpp[428];

    auto g_x_0_z_0_y_zz_z_x = buffer_1010_pdpp[429];

    auto g_x_0_z_0_y_zz_z_y = buffer_1010_pdpp[430];

    auto g_x_0_z_0_y_zz_z_z = buffer_1010_pdpp[431];

    auto g_x_0_z_0_z_xx_x_x = buffer_1010_pdpp[432];

    auto g_x_0_z_0_z_xx_x_y = buffer_1010_pdpp[433];

    auto g_x_0_z_0_z_xx_x_z = buffer_1010_pdpp[434];

    auto g_x_0_z_0_z_xx_y_x = buffer_1010_pdpp[435];

    auto g_x_0_z_0_z_xx_y_y = buffer_1010_pdpp[436];

    auto g_x_0_z_0_z_xx_y_z = buffer_1010_pdpp[437];

    auto g_x_0_z_0_z_xx_z_x = buffer_1010_pdpp[438];

    auto g_x_0_z_0_z_xx_z_y = buffer_1010_pdpp[439];

    auto g_x_0_z_0_z_xx_z_z = buffer_1010_pdpp[440];

    auto g_x_0_z_0_z_xy_x_x = buffer_1010_pdpp[441];

    auto g_x_0_z_0_z_xy_x_y = buffer_1010_pdpp[442];

    auto g_x_0_z_0_z_xy_x_z = buffer_1010_pdpp[443];

    auto g_x_0_z_0_z_xy_y_x = buffer_1010_pdpp[444];

    auto g_x_0_z_0_z_xy_y_y = buffer_1010_pdpp[445];

    auto g_x_0_z_0_z_xy_y_z = buffer_1010_pdpp[446];

    auto g_x_0_z_0_z_xy_z_x = buffer_1010_pdpp[447];

    auto g_x_0_z_0_z_xy_z_y = buffer_1010_pdpp[448];

    auto g_x_0_z_0_z_xy_z_z = buffer_1010_pdpp[449];

    auto g_x_0_z_0_z_xz_x_x = buffer_1010_pdpp[450];

    auto g_x_0_z_0_z_xz_x_y = buffer_1010_pdpp[451];

    auto g_x_0_z_0_z_xz_x_z = buffer_1010_pdpp[452];

    auto g_x_0_z_0_z_xz_y_x = buffer_1010_pdpp[453];

    auto g_x_0_z_0_z_xz_y_y = buffer_1010_pdpp[454];

    auto g_x_0_z_0_z_xz_y_z = buffer_1010_pdpp[455];

    auto g_x_0_z_0_z_xz_z_x = buffer_1010_pdpp[456];

    auto g_x_0_z_0_z_xz_z_y = buffer_1010_pdpp[457];

    auto g_x_0_z_0_z_xz_z_z = buffer_1010_pdpp[458];

    auto g_x_0_z_0_z_yy_x_x = buffer_1010_pdpp[459];

    auto g_x_0_z_0_z_yy_x_y = buffer_1010_pdpp[460];

    auto g_x_0_z_0_z_yy_x_z = buffer_1010_pdpp[461];

    auto g_x_0_z_0_z_yy_y_x = buffer_1010_pdpp[462];

    auto g_x_0_z_0_z_yy_y_y = buffer_1010_pdpp[463];

    auto g_x_0_z_0_z_yy_y_z = buffer_1010_pdpp[464];

    auto g_x_0_z_0_z_yy_z_x = buffer_1010_pdpp[465];

    auto g_x_0_z_0_z_yy_z_y = buffer_1010_pdpp[466];

    auto g_x_0_z_0_z_yy_z_z = buffer_1010_pdpp[467];

    auto g_x_0_z_0_z_yz_x_x = buffer_1010_pdpp[468];

    auto g_x_0_z_0_z_yz_x_y = buffer_1010_pdpp[469];

    auto g_x_0_z_0_z_yz_x_z = buffer_1010_pdpp[470];

    auto g_x_0_z_0_z_yz_y_x = buffer_1010_pdpp[471];

    auto g_x_0_z_0_z_yz_y_y = buffer_1010_pdpp[472];

    auto g_x_0_z_0_z_yz_y_z = buffer_1010_pdpp[473];

    auto g_x_0_z_0_z_yz_z_x = buffer_1010_pdpp[474];

    auto g_x_0_z_0_z_yz_z_y = buffer_1010_pdpp[475];

    auto g_x_0_z_0_z_yz_z_z = buffer_1010_pdpp[476];

    auto g_x_0_z_0_z_zz_x_x = buffer_1010_pdpp[477];

    auto g_x_0_z_0_z_zz_x_y = buffer_1010_pdpp[478];

    auto g_x_0_z_0_z_zz_x_z = buffer_1010_pdpp[479];

    auto g_x_0_z_0_z_zz_y_x = buffer_1010_pdpp[480];

    auto g_x_0_z_0_z_zz_y_y = buffer_1010_pdpp[481];

    auto g_x_0_z_0_z_zz_y_z = buffer_1010_pdpp[482];

    auto g_x_0_z_0_z_zz_z_x = buffer_1010_pdpp[483];

    auto g_x_0_z_0_z_zz_z_y = buffer_1010_pdpp[484];

    auto g_x_0_z_0_z_zz_z_z = buffer_1010_pdpp[485];

    auto g_y_0_x_0_x_xx_x_x = buffer_1010_pdpp[486];

    auto g_y_0_x_0_x_xx_x_y = buffer_1010_pdpp[487];

    auto g_y_0_x_0_x_xx_x_z = buffer_1010_pdpp[488];

    auto g_y_0_x_0_x_xx_y_x = buffer_1010_pdpp[489];

    auto g_y_0_x_0_x_xx_y_y = buffer_1010_pdpp[490];

    auto g_y_0_x_0_x_xx_y_z = buffer_1010_pdpp[491];

    auto g_y_0_x_0_x_xx_z_x = buffer_1010_pdpp[492];

    auto g_y_0_x_0_x_xx_z_y = buffer_1010_pdpp[493];

    auto g_y_0_x_0_x_xx_z_z = buffer_1010_pdpp[494];

    auto g_y_0_x_0_x_xy_x_x = buffer_1010_pdpp[495];

    auto g_y_0_x_0_x_xy_x_y = buffer_1010_pdpp[496];

    auto g_y_0_x_0_x_xy_x_z = buffer_1010_pdpp[497];

    auto g_y_0_x_0_x_xy_y_x = buffer_1010_pdpp[498];

    auto g_y_0_x_0_x_xy_y_y = buffer_1010_pdpp[499];

    auto g_y_0_x_0_x_xy_y_z = buffer_1010_pdpp[500];

    auto g_y_0_x_0_x_xy_z_x = buffer_1010_pdpp[501];

    auto g_y_0_x_0_x_xy_z_y = buffer_1010_pdpp[502];

    auto g_y_0_x_0_x_xy_z_z = buffer_1010_pdpp[503];

    auto g_y_0_x_0_x_xz_x_x = buffer_1010_pdpp[504];

    auto g_y_0_x_0_x_xz_x_y = buffer_1010_pdpp[505];

    auto g_y_0_x_0_x_xz_x_z = buffer_1010_pdpp[506];

    auto g_y_0_x_0_x_xz_y_x = buffer_1010_pdpp[507];

    auto g_y_0_x_0_x_xz_y_y = buffer_1010_pdpp[508];

    auto g_y_0_x_0_x_xz_y_z = buffer_1010_pdpp[509];

    auto g_y_0_x_0_x_xz_z_x = buffer_1010_pdpp[510];

    auto g_y_0_x_0_x_xz_z_y = buffer_1010_pdpp[511];

    auto g_y_0_x_0_x_xz_z_z = buffer_1010_pdpp[512];

    auto g_y_0_x_0_x_yy_x_x = buffer_1010_pdpp[513];

    auto g_y_0_x_0_x_yy_x_y = buffer_1010_pdpp[514];

    auto g_y_0_x_0_x_yy_x_z = buffer_1010_pdpp[515];

    auto g_y_0_x_0_x_yy_y_x = buffer_1010_pdpp[516];

    auto g_y_0_x_0_x_yy_y_y = buffer_1010_pdpp[517];

    auto g_y_0_x_0_x_yy_y_z = buffer_1010_pdpp[518];

    auto g_y_0_x_0_x_yy_z_x = buffer_1010_pdpp[519];

    auto g_y_0_x_0_x_yy_z_y = buffer_1010_pdpp[520];

    auto g_y_0_x_0_x_yy_z_z = buffer_1010_pdpp[521];

    auto g_y_0_x_0_x_yz_x_x = buffer_1010_pdpp[522];

    auto g_y_0_x_0_x_yz_x_y = buffer_1010_pdpp[523];

    auto g_y_0_x_0_x_yz_x_z = buffer_1010_pdpp[524];

    auto g_y_0_x_0_x_yz_y_x = buffer_1010_pdpp[525];

    auto g_y_0_x_0_x_yz_y_y = buffer_1010_pdpp[526];

    auto g_y_0_x_0_x_yz_y_z = buffer_1010_pdpp[527];

    auto g_y_0_x_0_x_yz_z_x = buffer_1010_pdpp[528];

    auto g_y_0_x_0_x_yz_z_y = buffer_1010_pdpp[529];

    auto g_y_0_x_0_x_yz_z_z = buffer_1010_pdpp[530];

    auto g_y_0_x_0_x_zz_x_x = buffer_1010_pdpp[531];

    auto g_y_0_x_0_x_zz_x_y = buffer_1010_pdpp[532];

    auto g_y_0_x_0_x_zz_x_z = buffer_1010_pdpp[533];

    auto g_y_0_x_0_x_zz_y_x = buffer_1010_pdpp[534];

    auto g_y_0_x_0_x_zz_y_y = buffer_1010_pdpp[535];

    auto g_y_0_x_0_x_zz_y_z = buffer_1010_pdpp[536];

    auto g_y_0_x_0_x_zz_z_x = buffer_1010_pdpp[537];

    auto g_y_0_x_0_x_zz_z_y = buffer_1010_pdpp[538];

    auto g_y_0_x_0_x_zz_z_z = buffer_1010_pdpp[539];

    auto g_y_0_x_0_y_xx_x_x = buffer_1010_pdpp[540];

    auto g_y_0_x_0_y_xx_x_y = buffer_1010_pdpp[541];

    auto g_y_0_x_0_y_xx_x_z = buffer_1010_pdpp[542];

    auto g_y_0_x_0_y_xx_y_x = buffer_1010_pdpp[543];

    auto g_y_0_x_0_y_xx_y_y = buffer_1010_pdpp[544];

    auto g_y_0_x_0_y_xx_y_z = buffer_1010_pdpp[545];

    auto g_y_0_x_0_y_xx_z_x = buffer_1010_pdpp[546];

    auto g_y_0_x_0_y_xx_z_y = buffer_1010_pdpp[547];

    auto g_y_0_x_0_y_xx_z_z = buffer_1010_pdpp[548];

    auto g_y_0_x_0_y_xy_x_x = buffer_1010_pdpp[549];

    auto g_y_0_x_0_y_xy_x_y = buffer_1010_pdpp[550];

    auto g_y_0_x_0_y_xy_x_z = buffer_1010_pdpp[551];

    auto g_y_0_x_0_y_xy_y_x = buffer_1010_pdpp[552];

    auto g_y_0_x_0_y_xy_y_y = buffer_1010_pdpp[553];

    auto g_y_0_x_0_y_xy_y_z = buffer_1010_pdpp[554];

    auto g_y_0_x_0_y_xy_z_x = buffer_1010_pdpp[555];

    auto g_y_0_x_0_y_xy_z_y = buffer_1010_pdpp[556];

    auto g_y_0_x_0_y_xy_z_z = buffer_1010_pdpp[557];

    auto g_y_0_x_0_y_xz_x_x = buffer_1010_pdpp[558];

    auto g_y_0_x_0_y_xz_x_y = buffer_1010_pdpp[559];

    auto g_y_0_x_0_y_xz_x_z = buffer_1010_pdpp[560];

    auto g_y_0_x_0_y_xz_y_x = buffer_1010_pdpp[561];

    auto g_y_0_x_0_y_xz_y_y = buffer_1010_pdpp[562];

    auto g_y_0_x_0_y_xz_y_z = buffer_1010_pdpp[563];

    auto g_y_0_x_0_y_xz_z_x = buffer_1010_pdpp[564];

    auto g_y_0_x_0_y_xz_z_y = buffer_1010_pdpp[565];

    auto g_y_0_x_0_y_xz_z_z = buffer_1010_pdpp[566];

    auto g_y_0_x_0_y_yy_x_x = buffer_1010_pdpp[567];

    auto g_y_0_x_0_y_yy_x_y = buffer_1010_pdpp[568];

    auto g_y_0_x_0_y_yy_x_z = buffer_1010_pdpp[569];

    auto g_y_0_x_0_y_yy_y_x = buffer_1010_pdpp[570];

    auto g_y_0_x_0_y_yy_y_y = buffer_1010_pdpp[571];

    auto g_y_0_x_0_y_yy_y_z = buffer_1010_pdpp[572];

    auto g_y_0_x_0_y_yy_z_x = buffer_1010_pdpp[573];

    auto g_y_0_x_0_y_yy_z_y = buffer_1010_pdpp[574];

    auto g_y_0_x_0_y_yy_z_z = buffer_1010_pdpp[575];

    auto g_y_0_x_0_y_yz_x_x = buffer_1010_pdpp[576];

    auto g_y_0_x_0_y_yz_x_y = buffer_1010_pdpp[577];

    auto g_y_0_x_0_y_yz_x_z = buffer_1010_pdpp[578];

    auto g_y_0_x_0_y_yz_y_x = buffer_1010_pdpp[579];

    auto g_y_0_x_0_y_yz_y_y = buffer_1010_pdpp[580];

    auto g_y_0_x_0_y_yz_y_z = buffer_1010_pdpp[581];

    auto g_y_0_x_0_y_yz_z_x = buffer_1010_pdpp[582];

    auto g_y_0_x_0_y_yz_z_y = buffer_1010_pdpp[583];

    auto g_y_0_x_0_y_yz_z_z = buffer_1010_pdpp[584];

    auto g_y_0_x_0_y_zz_x_x = buffer_1010_pdpp[585];

    auto g_y_0_x_0_y_zz_x_y = buffer_1010_pdpp[586];

    auto g_y_0_x_0_y_zz_x_z = buffer_1010_pdpp[587];

    auto g_y_0_x_0_y_zz_y_x = buffer_1010_pdpp[588];

    auto g_y_0_x_0_y_zz_y_y = buffer_1010_pdpp[589];

    auto g_y_0_x_0_y_zz_y_z = buffer_1010_pdpp[590];

    auto g_y_0_x_0_y_zz_z_x = buffer_1010_pdpp[591];

    auto g_y_0_x_0_y_zz_z_y = buffer_1010_pdpp[592];

    auto g_y_0_x_0_y_zz_z_z = buffer_1010_pdpp[593];

    auto g_y_0_x_0_z_xx_x_x = buffer_1010_pdpp[594];

    auto g_y_0_x_0_z_xx_x_y = buffer_1010_pdpp[595];

    auto g_y_0_x_0_z_xx_x_z = buffer_1010_pdpp[596];

    auto g_y_0_x_0_z_xx_y_x = buffer_1010_pdpp[597];

    auto g_y_0_x_0_z_xx_y_y = buffer_1010_pdpp[598];

    auto g_y_0_x_0_z_xx_y_z = buffer_1010_pdpp[599];

    auto g_y_0_x_0_z_xx_z_x = buffer_1010_pdpp[600];

    auto g_y_0_x_0_z_xx_z_y = buffer_1010_pdpp[601];

    auto g_y_0_x_0_z_xx_z_z = buffer_1010_pdpp[602];

    auto g_y_0_x_0_z_xy_x_x = buffer_1010_pdpp[603];

    auto g_y_0_x_0_z_xy_x_y = buffer_1010_pdpp[604];

    auto g_y_0_x_0_z_xy_x_z = buffer_1010_pdpp[605];

    auto g_y_0_x_0_z_xy_y_x = buffer_1010_pdpp[606];

    auto g_y_0_x_0_z_xy_y_y = buffer_1010_pdpp[607];

    auto g_y_0_x_0_z_xy_y_z = buffer_1010_pdpp[608];

    auto g_y_0_x_0_z_xy_z_x = buffer_1010_pdpp[609];

    auto g_y_0_x_0_z_xy_z_y = buffer_1010_pdpp[610];

    auto g_y_0_x_0_z_xy_z_z = buffer_1010_pdpp[611];

    auto g_y_0_x_0_z_xz_x_x = buffer_1010_pdpp[612];

    auto g_y_0_x_0_z_xz_x_y = buffer_1010_pdpp[613];

    auto g_y_0_x_0_z_xz_x_z = buffer_1010_pdpp[614];

    auto g_y_0_x_0_z_xz_y_x = buffer_1010_pdpp[615];

    auto g_y_0_x_0_z_xz_y_y = buffer_1010_pdpp[616];

    auto g_y_0_x_0_z_xz_y_z = buffer_1010_pdpp[617];

    auto g_y_0_x_0_z_xz_z_x = buffer_1010_pdpp[618];

    auto g_y_0_x_0_z_xz_z_y = buffer_1010_pdpp[619];

    auto g_y_0_x_0_z_xz_z_z = buffer_1010_pdpp[620];

    auto g_y_0_x_0_z_yy_x_x = buffer_1010_pdpp[621];

    auto g_y_0_x_0_z_yy_x_y = buffer_1010_pdpp[622];

    auto g_y_0_x_0_z_yy_x_z = buffer_1010_pdpp[623];

    auto g_y_0_x_0_z_yy_y_x = buffer_1010_pdpp[624];

    auto g_y_0_x_0_z_yy_y_y = buffer_1010_pdpp[625];

    auto g_y_0_x_0_z_yy_y_z = buffer_1010_pdpp[626];

    auto g_y_0_x_0_z_yy_z_x = buffer_1010_pdpp[627];

    auto g_y_0_x_0_z_yy_z_y = buffer_1010_pdpp[628];

    auto g_y_0_x_0_z_yy_z_z = buffer_1010_pdpp[629];

    auto g_y_0_x_0_z_yz_x_x = buffer_1010_pdpp[630];

    auto g_y_0_x_0_z_yz_x_y = buffer_1010_pdpp[631];

    auto g_y_0_x_0_z_yz_x_z = buffer_1010_pdpp[632];

    auto g_y_0_x_0_z_yz_y_x = buffer_1010_pdpp[633];

    auto g_y_0_x_0_z_yz_y_y = buffer_1010_pdpp[634];

    auto g_y_0_x_0_z_yz_y_z = buffer_1010_pdpp[635];

    auto g_y_0_x_0_z_yz_z_x = buffer_1010_pdpp[636];

    auto g_y_0_x_0_z_yz_z_y = buffer_1010_pdpp[637];

    auto g_y_0_x_0_z_yz_z_z = buffer_1010_pdpp[638];

    auto g_y_0_x_0_z_zz_x_x = buffer_1010_pdpp[639];

    auto g_y_0_x_0_z_zz_x_y = buffer_1010_pdpp[640];

    auto g_y_0_x_0_z_zz_x_z = buffer_1010_pdpp[641];

    auto g_y_0_x_0_z_zz_y_x = buffer_1010_pdpp[642];

    auto g_y_0_x_0_z_zz_y_y = buffer_1010_pdpp[643];

    auto g_y_0_x_0_z_zz_y_z = buffer_1010_pdpp[644];

    auto g_y_0_x_0_z_zz_z_x = buffer_1010_pdpp[645];

    auto g_y_0_x_0_z_zz_z_y = buffer_1010_pdpp[646];

    auto g_y_0_x_0_z_zz_z_z = buffer_1010_pdpp[647];

    auto g_y_0_y_0_x_xx_x_x = buffer_1010_pdpp[648];

    auto g_y_0_y_0_x_xx_x_y = buffer_1010_pdpp[649];

    auto g_y_0_y_0_x_xx_x_z = buffer_1010_pdpp[650];

    auto g_y_0_y_0_x_xx_y_x = buffer_1010_pdpp[651];

    auto g_y_0_y_0_x_xx_y_y = buffer_1010_pdpp[652];

    auto g_y_0_y_0_x_xx_y_z = buffer_1010_pdpp[653];

    auto g_y_0_y_0_x_xx_z_x = buffer_1010_pdpp[654];

    auto g_y_0_y_0_x_xx_z_y = buffer_1010_pdpp[655];

    auto g_y_0_y_0_x_xx_z_z = buffer_1010_pdpp[656];

    auto g_y_0_y_0_x_xy_x_x = buffer_1010_pdpp[657];

    auto g_y_0_y_0_x_xy_x_y = buffer_1010_pdpp[658];

    auto g_y_0_y_0_x_xy_x_z = buffer_1010_pdpp[659];

    auto g_y_0_y_0_x_xy_y_x = buffer_1010_pdpp[660];

    auto g_y_0_y_0_x_xy_y_y = buffer_1010_pdpp[661];

    auto g_y_0_y_0_x_xy_y_z = buffer_1010_pdpp[662];

    auto g_y_0_y_0_x_xy_z_x = buffer_1010_pdpp[663];

    auto g_y_0_y_0_x_xy_z_y = buffer_1010_pdpp[664];

    auto g_y_0_y_0_x_xy_z_z = buffer_1010_pdpp[665];

    auto g_y_0_y_0_x_xz_x_x = buffer_1010_pdpp[666];

    auto g_y_0_y_0_x_xz_x_y = buffer_1010_pdpp[667];

    auto g_y_0_y_0_x_xz_x_z = buffer_1010_pdpp[668];

    auto g_y_0_y_0_x_xz_y_x = buffer_1010_pdpp[669];

    auto g_y_0_y_0_x_xz_y_y = buffer_1010_pdpp[670];

    auto g_y_0_y_0_x_xz_y_z = buffer_1010_pdpp[671];

    auto g_y_0_y_0_x_xz_z_x = buffer_1010_pdpp[672];

    auto g_y_0_y_0_x_xz_z_y = buffer_1010_pdpp[673];

    auto g_y_0_y_0_x_xz_z_z = buffer_1010_pdpp[674];

    auto g_y_0_y_0_x_yy_x_x = buffer_1010_pdpp[675];

    auto g_y_0_y_0_x_yy_x_y = buffer_1010_pdpp[676];

    auto g_y_0_y_0_x_yy_x_z = buffer_1010_pdpp[677];

    auto g_y_0_y_0_x_yy_y_x = buffer_1010_pdpp[678];

    auto g_y_0_y_0_x_yy_y_y = buffer_1010_pdpp[679];

    auto g_y_0_y_0_x_yy_y_z = buffer_1010_pdpp[680];

    auto g_y_0_y_0_x_yy_z_x = buffer_1010_pdpp[681];

    auto g_y_0_y_0_x_yy_z_y = buffer_1010_pdpp[682];

    auto g_y_0_y_0_x_yy_z_z = buffer_1010_pdpp[683];

    auto g_y_0_y_0_x_yz_x_x = buffer_1010_pdpp[684];

    auto g_y_0_y_0_x_yz_x_y = buffer_1010_pdpp[685];

    auto g_y_0_y_0_x_yz_x_z = buffer_1010_pdpp[686];

    auto g_y_0_y_0_x_yz_y_x = buffer_1010_pdpp[687];

    auto g_y_0_y_0_x_yz_y_y = buffer_1010_pdpp[688];

    auto g_y_0_y_0_x_yz_y_z = buffer_1010_pdpp[689];

    auto g_y_0_y_0_x_yz_z_x = buffer_1010_pdpp[690];

    auto g_y_0_y_0_x_yz_z_y = buffer_1010_pdpp[691];

    auto g_y_0_y_0_x_yz_z_z = buffer_1010_pdpp[692];

    auto g_y_0_y_0_x_zz_x_x = buffer_1010_pdpp[693];

    auto g_y_0_y_0_x_zz_x_y = buffer_1010_pdpp[694];

    auto g_y_0_y_0_x_zz_x_z = buffer_1010_pdpp[695];

    auto g_y_0_y_0_x_zz_y_x = buffer_1010_pdpp[696];

    auto g_y_0_y_0_x_zz_y_y = buffer_1010_pdpp[697];

    auto g_y_0_y_0_x_zz_y_z = buffer_1010_pdpp[698];

    auto g_y_0_y_0_x_zz_z_x = buffer_1010_pdpp[699];

    auto g_y_0_y_0_x_zz_z_y = buffer_1010_pdpp[700];

    auto g_y_0_y_0_x_zz_z_z = buffer_1010_pdpp[701];

    auto g_y_0_y_0_y_xx_x_x = buffer_1010_pdpp[702];

    auto g_y_0_y_0_y_xx_x_y = buffer_1010_pdpp[703];

    auto g_y_0_y_0_y_xx_x_z = buffer_1010_pdpp[704];

    auto g_y_0_y_0_y_xx_y_x = buffer_1010_pdpp[705];

    auto g_y_0_y_0_y_xx_y_y = buffer_1010_pdpp[706];

    auto g_y_0_y_0_y_xx_y_z = buffer_1010_pdpp[707];

    auto g_y_0_y_0_y_xx_z_x = buffer_1010_pdpp[708];

    auto g_y_0_y_0_y_xx_z_y = buffer_1010_pdpp[709];

    auto g_y_0_y_0_y_xx_z_z = buffer_1010_pdpp[710];

    auto g_y_0_y_0_y_xy_x_x = buffer_1010_pdpp[711];

    auto g_y_0_y_0_y_xy_x_y = buffer_1010_pdpp[712];

    auto g_y_0_y_0_y_xy_x_z = buffer_1010_pdpp[713];

    auto g_y_0_y_0_y_xy_y_x = buffer_1010_pdpp[714];

    auto g_y_0_y_0_y_xy_y_y = buffer_1010_pdpp[715];

    auto g_y_0_y_0_y_xy_y_z = buffer_1010_pdpp[716];

    auto g_y_0_y_0_y_xy_z_x = buffer_1010_pdpp[717];

    auto g_y_0_y_0_y_xy_z_y = buffer_1010_pdpp[718];

    auto g_y_0_y_0_y_xy_z_z = buffer_1010_pdpp[719];

    auto g_y_0_y_0_y_xz_x_x = buffer_1010_pdpp[720];

    auto g_y_0_y_0_y_xz_x_y = buffer_1010_pdpp[721];

    auto g_y_0_y_0_y_xz_x_z = buffer_1010_pdpp[722];

    auto g_y_0_y_0_y_xz_y_x = buffer_1010_pdpp[723];

    auto g_y_0_y_0_y_xz_y_y = buffer_1010_pdpp[724];

    auto g_y_0_y_0_y_xz_y_z = buffer_1010_pdpp[725];

    auto g_y_0_y_0_y_xz_z_x = buffer_1010_pdpp[726];

    auto g_y_0_y_0_y_xz_z_y = buffer_1010_pdpp[727];

    auto g_y_0_y_0_y_xz_z_z = buffer_1010_pdpp[728];

    auto g_y_0_y_0_y_yy_x_x = buffer_1010_pdpp[729];

    auto g_y_0_y_0_y_yy_x_y = buffer_1010_pdpp[730];

    auto g_y_0_y_0_y_yy_x_z = buffer_1010_pdpp[731];

    auto g_y_0_y_0_y_yy_y_x = buffer_1010_pdpp[732];

    auto g_y_0_y_0_y_yy_y_y = buffer_1010_pdpp[733];

    auto g_y_0_y_0_y_yy_y_z = buffer_1010_pdpp[734];

    auto g_y_0_y_0_y_yy_z_x = buffer_1010_pdpp[735];

    auto g_y_0_y_0_y_yy_z_y = buffer_1010_pdpp[736];

    auto g_y_0_y_0_y_yy_z_z = buffer_1010_pdpp[737];

    auto g_y_0_y_0_y_yz_x_x = buffer_1010_pdpp[738];

    auto g_y_0_y_0_y_yz_x_y = buffer_1010_pdpp[739];

    auto g_y_0_y_0_y_yz_x_z = buffer_1010_pdpp[740];

    auto g_y_0_y_0_y_yz_y_x = buffer_1010_pdpp[741];

    auto g_y_0_y_0_y_yz_y_y = buffer_1010_pdpp[742];

    auto g_y_0_y_0_y_yz_y_z = buffer_1010_pdpp[743];

    auto g_y_0_y_0_y_yz_z_x = buffer_1010_pdpp[744];

    auto g_y_0_y_0_y_yz_z_y = buffer_1010_pdpp[745];

    auto g_y_0_y_0_y_yz_z_z = buffer_1010_pdpp[746];

    auto g_y_0_y_0_y_zz_x_x = buffer_1010_pdpp[747];

    auto g_y_0_y_0_y_zz_x_y = buffer_1010_pdpp[748];

    auto g_y_0_y_0_y_zz_x_z = buffer_1010_pdpp[749];

    auto g_y_0_y_0_y_zz_y_x = buffer_1010_pdpp[750];

    auto g_y_0_y_0_y_zz_y_y = buffer_1010_pdpp[751];

    auto g_y_0_y_0_y_zz_y_z = buffer_1010_pdpp[752];

    auto g_y_0_y_0_y_zz_z_x = buffer_1010_pdpp[753];

    auto g_y_0_y_0_y_zz_z_y = buffer_1010_pdpp[754];

    auto g_y_0_y_0_y_zz_z_z = buffer_1010_pdpp[755];

    auto g_y_0_y_0_z_xx_x_x = buffer_1010_pdpp[756];

    auto g_y_0_y_0_z_xx_x_y = buffer_1010_pdpp[757];

    auto g_y_0_y_0_z_xx_x_z = buffer_1010_pdpp[758];

    auto g_y_0_y_0_z_xx_y_x = buffer_1010_pdpp[759];

    auto g_y_0_y_0_z_xx_y_y = buffer_1010_pdpp[760];

    auto g_y_0_y_0_z_xx_y_z = buffer_1010_pdpp[761];

    auto g_y_0_y_0_z_xx_z_x = buffer_1010_pdpp[762];

    auto g_y_0_y_0_z_xx_z_y = buffer_1010_pdpp[763];

    auto g_y_0_y_0_z_xx_z_z = buffer_1010_pdpp[764];

    auto g_y_0_y_0_z_xy_x_x = buffer_1010_pdpp[765];

    auto g_y_0_y_0_z_xy_x_y = buffer_1010_pdpp[766];

    auto g_y_0_y_0_z_xy_x_z = buffer_1010_pdpp[767];

    auto g_y_0_y_0_z_xy_y_x = buffer_1010_pdpp[768];

    auto g_y_0_y_0_z_xy_y_y = buffer_1010_pdpp[769];

    auto g_y_0_y_0_z_xy_y_z = buffer_1010_pdpp[770];

    auto g_y_0_y_0_z_xy_z_x = buffer_1010_pdpp[771];

    auto g_y_0_y_0_z_xy_z_y = buffer_1010_pdpp[772];

    auto g_y_0_y_0_z_xy_z_z = buffer_1010_pdpp[773];

    auto g_y_0_y_0_z_xz_x_x = buffer_1010_pdpp[774];

    auto g_y_0_y_0_z_xz_x_y = buffer_1010_pdpp[775];

    auto g_y_0_y_0_z_xz_x_z = buffer_1010_pdpp[776];

    auto g_y_0_y_0_z_xz_y_x = buffer_1010_pdpp[777];

    auto g_y_0_y_0_z_xz_y_y = buffer_1010_pdpp[778];

    auto g_y_0_y_0_z_xz_y_z = buffer_1010_pdpp[779];

    auto g_y_0_y_0_z_xz_z_x = buffer_1010_pdpp[780];

    auto g_y_0_y_0_z_xz_z_y = buffer_1010_pdpp[781];

    auto g_y_0_y_0_z_xz_z_z = buffer_1010_pdpp[782];

    auto g_y_0_y_0_z_yy_x_x = buffer_1010_pdpp[783];

    auto g_y_0_y_0_z_yy_x_y = buffer_1010_pdpp[784];

    auto g_y_0_y_0_z_yy_x_z = buffer_1010_pdpp[785];

    auto g_y_0_y_0_z_yy_y_x = buffer_1010_pdpp[786];

    auto g_y_0_y_0_z_yy_y_y = buffer_1010_pdpp[787];

    auto g_y_0_y_0_z_yy_y_z = buffer_1010_pdpp[788];

    auto g_y_0_y_0_z_yy_z_x = buffer_1010_pdpp[789];

    auto g_y_0_y_0_z_yy_z_y = buffer_1010_pdpp[790];

    auto g_y_0_y_0_z_yy_z_z = buffer_1010_pdpp[791];

    auto g_y_0_y_0_z_yz_x_x = buffer_1010_pdpp[792];

    auto g_y_0_y_0_z_yz_x_y = buffer_1010_pdpp[793];

    auto g_y_0_y_0_z_yz_x_z = buffer_1010_pdpp[794];

    auto g_y_0_y_0_z_yz_y_x = buffer_1010_pdpp[795];

    auto g_y_0_y_0_z_yz_y_y = buffer_1010_pdpp[796];

    auto g_y_0_y_0_z_yz_y_z = buffer_1010_pdpp[797];

    auto g_y_0_y_0_z_yz_z_x = buffer_1010_pdpp[798];

    auto g_y_0_y_0_z_yz_z_y = buffer_1010_pdpp[799];

    auto g_y_0_y_0_z_yz_z_z = buffer_1010_pdpp[800];

    auto g_y_0_y_0_z_zz_x_x = buffer_1010_pdpp[801];

    auto g_y_0_y_0_z_zz_x_y = buffer_1010_pdpp[802];

    auto g_y_0_y_0_z_zz_x_z = buffer_1010_pdpp[803];

    auto g_y_0_y_0_z_zz_y_x = buffer_1010_pdpp[804];

    auto g_y_0_y_0_z_zz_y_y = buffer_1010_pdpp[805];

    auto g_y_0_y_0_z_zz_y_z = buffer_1010_pdpp[806];

    auto g_y_0_y_0_z_zz_z_x = buffer_1010_pdpp[807];

    auto g_y_0_y_0_z_zz_z_y = buffer_1010_pdpp[808];

    auto g_y_0_y_0_z_zz_z_z = buffer_1010_pdpp[809];

    auto g_y_0_z_0_x_xx_x_x = buffer_1010_pdpp[810];

    auto g_y_0_z_0_x_xx_x_y = buffer_1010_pdpp[811];

    auto g_y_0_z_0_x_xx_x_z = buffer_1010_pdpp[812];

    auto g_y_0_z_0_x_xx_y_x = buffer_1010_pdpp[813];

    auto g_y_0_z_0_x_xx_y_y = buffer_1010_pdpp[814];

    auto g_y_0_z_0_x_xx_y_z = buffer_1010_pdpp[815];

    auto g_y_0_z_0_x_xx_z_x = buffer_1010_pdpp[816];

    auto g_y_0_z_0_x_xx_z_y = buffer_1010_pdpp[817];

    auto g_y_0_z_0_x_xx_z_z = buffer_1010_pdpp[818];

    auto g_y_0_z_0_x_xy_x_x = buffer_1010_pdpp[819];

    auto g_y_0_z_0_x_xy_x_y = buffer_1010_pdpp[820];

    auto g_y_0_z_0_x_xy_x_z = buffer_1010_pdpp[821];

    auto g_y_0_z_0_x_xy_y_x = buffer_1010_pdpp[822];

    auto g_y_0_z_0_x_xy_y_y = buffer_1010_pdpp[823];

    auto g_y_0_z_0_x_xy_y_z = buffer_1010_pdpp[824];

    auto g_y_0_z_0_x_xy_z_x = buffer_1010_pdpp[825];

    auto g_y_0_z_0_x_xy_z_y = buffer_1010_pdpp[826];

    auto g_y_0_z_0_x_xy_z_z = buffer_1010_pdpp[827];

    auto g_y_0_z_0_x_xz_x_x = buffer_1010_pdpp[828];

    auto g_y_0_z_0_x_xz_x_y = buffer_1010_pdpp[829];

    auto g_y_0_z_0_x_xz_x_z = buffer_1010_pdpp[830];

    auto g_y_0_z_0_x_xz_y_x = buffer_1010_pdpp[831];

    auto g_y_0_z_0_x_xz_y_y = buffer_1010_pdpp[832];

    auto g_y_0_z_0_x_xz_y_z = buffer_1010_pdpp[833];

    auto g_y_0_z_0_x_xz_z_x = buffer_1010_pdpp[834];

    auto g_y_0_z_0_x_xz_z_y = buffer_1010_pdpp[835];

    auto g_y_0_z_0_x_xz_z_z = buffer_1010_pdpp[836];

    auto g_y_0_z_0_x_yy_x_x = buffer_1010_pdpp[837];

    auto g_y_0_z_0_x_yy_x_y = buffer_1010_pdpp[838];

    auto g_y_0_z_0_x_yy_x_z = buffer_1010_pdpp[839];

    auto g_y_0_z_0_x_yy_y_x = buffer_1010_pdpp[840];

    auto g_y_0_z_0_x_yy_y_y = buffer_1010_pdpp[841];

    auto g_y_0_z_0_x_yy_y_z = buffer_1010_pdpp[842];

    auto g_y_0_z_0_x_yy_z_x = buffer_1010_pdpp[843];

    auto g_y_0_z_0_x_yy_z_y = buffer_1010_pdpp[844];

    auto g_y_0_z_0_x_yy_z_z = buffer_1010_pdpp[845];

    auto g_y_0_z_0_x_yz_x_x = buffer_1010_pdpp[846];

    auto g_y_0_z_0_x_yz_x_y = buffer_1010_pdpp[847];

    auto g_y_0_z_0_x_yz_x_z = buffer_1010_pdpp[848];

    auto g_y_0_z_0_x_yz_y_x = buffer_1010_pdpp[849];

    auto g_y_0_z_0_x_yz_y_y = buffer_1010_pdpp[850];

    auto g_y_0_z_0_x_yz_y_z = buffer_1010_pdpp[851];

    auto g_y_0_z_0_x_yz_z_x = buffer_1010_pdpp[852];

    auto g_y_0_z_0_x_yz_z_y = buffer_1010_pdpp[853];

    auto g_y_0_z_0_x_yz_z_z = buffer_1010_pdpp[854];

    auto g_y_0_z_0_x_zz_x_x = buffer_1010_pdpp[855];

    auto g_y_0_z_0_x_zz_x_y = buffer_1010_pdpp[856];

    auto g_y_0_z_0_x_zz_x_z = buffer_1010_pdpp[857];

    auto g_y_0_z_0_x_zz_y_x = buffer_1010_pdpp[858];

    auto g_y_0_z_0_x_zz_y_y = buffer_1010_pdpp[859];

    auto g_y_0_z_0_x_zz_y_z = buffer_1010_pdpp[860];

    auto g_y_0_z_0_x_zz_z_x = buffer_1010_pdpp[861];

    auto g_y_0_z_0_x_zz_z_y = buffer_1010_pdpp[862];

    auto g_y_0_z_0_x_zz_z_z = buffer_1010_pdpp[863];

    auto g_y_0_z_0_y_xx_x_x = buffer_1010_pdpp[864];

    auto g_y_0_z_0_y_xx_x_y = buffer_1010_pdpp[865];

    auto g_y_0_z_0_y_xx_x_z = buffer_1010_pdpp[866];

    auto g_y_0_z_0_y_xx_y_x = buffer_1010_pdpp[867];

    auto g_y_0_z_0_y_xx_y_y = buffer_1010_pdpp[868];

    auto g_y_0_z_0_y_xx_y_z = buffer_1010_pdpp[869];

    auto g_y_0_z_0_y_xx_z_x = buffer_1010_pdpp[870];

    auto g_y_0_z_0_y_xx_z_y = buffer_1010_pdpp[871];

    auto g_y_0_z_0_y_xx_z_z = buffer_1010_pdpp[872];

    auto g_y_0_z_0_y_xy_x_x = buffer_1010_pdpp[873];

    auto g_y_0_z_0_y_xy_x_y = buffer_1010_pdpp[874];

    auto g_y_0_z_0_y_xy_x_z = buffer_1010_pdpp[875];

    auto g_y_0_z_0_y_xy_y_x = buffer_1010_pdpp[876];

    auto g_y_0_z_0_y_xy_y_y = buffer_1010_pdpp[877];

    auto g_y_0_z_0_y_xy_y_z = buffer_1010_pdpp[878];

    auto g_y_0_z_0_y_xy_z_x = buffer_1010_pdpp[879];

    auto g_y_0_z_0_y_xy_z_y = buffer_1010_pdpp[880];

    auto g_y_0_z_0_y_xy_z_z = buffer_1010_pdpp[881];

    auto g_y_0_z_0_y_xz_x_x = buffer_1010_pdpp[882];

    auto g_y_0_z_0_y_xz_x_y = buffer_1010_pdpp[883];

    auto g_y_0_z_0_y_xz_x_z = buffer_1010_pdpp[884];

    auto g_y_0_z_0_y_xz_y_x = buffer_1010_pdpp[885];

    auto g_y_0_z_0_y_xz_y_y = buffer_1010_pdpp[886];

    auto g_y_0_z_0_y_xz_y_z = buffer_1010_pdpp[887];

    auto g_y_0_z_0_y_xz_z_x = buffer_1010_pdpp[888];

    auto g_y_0_z_0_y_xz_z_y = buffer_1010_pdpp[889];

    auto g_y_0_z_0_y_xz_z_z = buffer_1010_pdpp[890];

    auto g_y_0_z_0_y_yy_x_x = buffer_1010_pdpp[891];

    auto g_y_0_z_0_y_yy_x_y = buffer_1010_pdpp[892];

    auto g_y_0_z_0_y_yy_x_z = buffer_1010_pdpp[893];

    auto g_y_0_z_0_y_yy_y_x = buffer_1010_pdpp[894];

    auto g_y_0_z_0_y_yy_y_y = buffer_1010_pdpp[895];

    auto g_y_0_z_0_y_yy_y_z = buffer_1010_pdpp[896];

    auto g_y_0_z_0_y_yy_z_x = buffer_1010_pdpp[897];

    auto g_y_0_z_0_y_yy_z_y = buffer_1010_pdpp[898];

    auto g_y_0_z_0_y_yy_z_z = buffer_1010_pdpp[899];

    auto g_y_0_z_0_y_yz_x_x = buffer_1010_pdpp[900];

    auto g_y_0_z_0_y_yz_x_y = buffer_1010_pdpp[901];

    auto g_y_0_z_0_y_yz_x_z = buffer_1010_pdpp[902];

    auto g_y_0_z_0_y_yz_y_x = buffer_1010_pdpp[903];

    auto g_y_0_z_0_y_yz_y_y = buffer_1010_pdpp[904];

    auto g_y_0_z_0_y_yz_y_z = buffer_1010_pdpp[905];

    auto g_y_0_z_0_y_yz_z_x = buffer_1010_pdpp[906];

    auto g_y_0_z_0_y_yz_z_y = buffer_1010_pdpp[907];

    auto g_y_0_z_0_y_yz_z_z = buffer_1010_pdpp[908];

    auto g_y_0_z_0_y_zz_x_x = buffer_1010_pdpp[909];

    auto g_y_0_z_0_y_zz_x_y = buffer_1010_pdpp[910];

    auto g_y_0_z_0_y_zz_x_z = buffer_1010_pdpp[911];

    auto g_y_0_z_0_y_zz_y_x = buffer_1010_pdpp[912];

    auto g_y_0_z_0_y_zz_y_y = buffer_1010_pdpp[913];

    auto g_y_0_z_0_y_zz_y_z = buffer_1010_pdpp[914];

    auto g_y_0_z_0_y_zz_z_x = buffer_1010_pdpp[915];

    auto g_y_0_z_0_y_zz_z_y = buffer_1010_pdpp[916];

    auto g_y_0_z_0_y_zz_z_z = buffer_1010_pdpp[917];

    auto g_y_0_z_0_z_xx_x_x = buffer_1010_pdpp[918];

    auto g_y_0_z_0_z_xx_x_y = buffer_1010_pdpp[919];

    auto g_y_0_z_0_z_xx_x_z = buffer_1010_pdpp[920];

    auto g_y_0_z_0_z_xx_y_x = buffer_1010_pdpp[921];

    auto g_y_0_z_0_z_xx_y_y = buffer_1010_pdpp[922];

    auto g_y_0_z_0_z_xx_y_z = buffer_1010_pdpp[923];

    auto g_y_0_z_0_z_xx_z_x = buffer_1010_pdpp[924];

    auto g_y_0_z_0_z_xx_z_y = buffer_1010_pdpp[925];

    auto g_y_0_z_0_z_xx_z_z = buffer_1010_pdpp[926];

    auto g_y_0_z_0_z_xy_x_x = buffer_1010_pdpp[927];

    auto g_y_0_z_0_z_xy_x_y = buffer_1010_pdpp[928];

    auto g_y_0_z_0_z_xy_x_z = buffer_1010_pdpp[929];

    auto g_y_0_z_0_z_xy_y_x = buffer_1010_pdpp[930];

    auto g_y_0_z_0_z_xy_y_y = buffer_1010_pdpp[931];

    auto g_y_0_z_0_z_xy_y_z = buffer_1010_pdpp[932];

    auto g_y_0_z_0_z_xy_z_x = buffer_1010_pdpp[933];

    auto g_y_0_z_0_z_xy_z_y = buffer_1010_pdpp[934];

    auto g_y_0_z_0_z_xy_z_z = buffer_1010_pdpp[935];

    auto g_y_0_z_0_z_xz_x_x = buffer_1010_pdpp[936];

    auto g_y_0_z_0_z_xz_x_y = buffer_1010_pdpp[937];

    auto g_y_0_z_0_z_xz_x_z = buffer_1010_pdpp[938];

    auto g_y_0_z_0_z_xz_y_x = buffer_1010_pdpp[939];

    auto g_y_0_z_0_z_xz_y_y = buffer_1010_pdpp[940];

    auto g_y_0_z_0_z_xz_y_z = buffer_1010_pdpp[941];

    auto g_y_0_z_0_z_xz_z_x = buffer_1010_pdpp[942];

    auto g_y_0_z_0_z_xz_z_y = buffer_1010_pdpp[943];

    auto g_y_0_z_0_z_xz_z_z = buffer_1010_pdpp[944];

    auto g_y_0_z_0_z_yy_x_x = buffer_1010_pdpp[945];

    auto g_y_0_z_0_z_yy_x_y = buffer_1010_pdpp[946];

    auto g_y_0_z_0_z_yy_x_z = buffer_1010_pdpp[947];

    auto g_y_0_z_0_z_yy_y_x = buffer_1010_pdpp[948];

    auto g_y_0_z_0_z_yy_y_y = buffer_1010_pdpp[949];

    auto g_y_0_z_0_z_yy_y_z = buffer_1010_pdpp[950];

    auto g_y_0_z_0_z_yy_z_x = buffer_1010_pdpp[951];

    auto g_y_0_z_0_z_yy_z_y = buffer_1010_pdpp[952];

    auto g_y_0_z_0_z_yy_z_z = buffer_1010_pdpp[953];

    auto g_y_0_z_0_z_yz_x_x = buffer_1010_pdpp[954];

    auto g_y_0_z_0_z_yz_x_y = buffer_1010_pdpp[955];

    auto g_y_0_z_0_z_yz_x_z = buffer_1010_pdpp[956];

    auto g_y_0_z_0_z_yz_y_x = buffer_1010_pdpp[957];

    auto g_y_0_z_0_z_yz_y_y = buffer_1010_pdpp[958];

    auto g_y_0_z_0_z_yz_y_z = buffer_1010_pdpp[959];

    auto g_y_0_z_0_z_yz_z_x = buffer_1010_pdpp[960];

    auto g_y_0_z_0_z_yz_z_y = buffer_1010_pdpp[961];

    auto g_y_0_z_0_z_yz_z_z = buffer_1010_pdpp[962];

    auto g_y_0_z_0_z_zz_x_x = buffer_1010_pdpp[963];

    auto g_y_0_z_0_z_zz_x_y = buffer_1010_pdpp[964];

    auto g_y_0_z_0_z_zz_x_z = buffer_1010_pdpp[965];

    auto g_y_0_z_0_z_zz_y_x = buffer_1010_pdpp[966];

    auto g_y_0_z_0_z_zz_y_y = buffer_1010_pdpp[967];

    auto g_y_0_z_0_z_zz_y_z = buffer_1010_pdpp[968];

    auto g_y_0_z_0_z_zz_z_x = buffer_1010_pdpp[969];

    auto g_y_0_z_0_z_zz_z_y = buffer_1010_pdpp[970];

    auto g_y_0_z_0_z_zz_z_z = buffer_1010_pdpp[971];

    auto g_z_0_x_0_x_xx_x_x = buffer_1010_pdpp[972];

    auto g_z_0_x_0_x_xx_x_y = buffer_1010_pdpp[973];

    auto g_z_0_x_0_x_xx_x_z = buffer_1010_pdpp[974];

    auto g_z_0_x_0_x_xx_y_x = buffer_1010_pdpp[975];

    auto g_z_0_x_0_x_xx_y_y = buffer_1010_pdpp[976];

    auto g_z_0_x_0_x_xx_y_z = buffer_1010_pdpp[977];

    auto g_z_0_x_0_x_xx_z_x = buffer_1010_pdpp[978];

    auto g_z_0_x_0_x_xx_z_y = buffer_1010_pdpp[979];

    auto g_z_0_x_0_x_xx_z_z = buffer_1010_pdpp[980];

    auto g_z_0_x_0_x_xy_x_x = buffer_1010_pdpp[981];

    auto g_z_0_x_0_x_xy_x_y = buffer_1010_pdpp[982];

    auto g_z_0_x_0_x_xy_x_z = buffer_1010_pdpp[983];

    auto g_z_0_x_0_x_xy_y_x = buffer_1010_pdpp[984];

    auto g_z_0_x_0_x_xy_y_y = buffer_1010_pdpp[985];

    auto g_z_0_x_0_x_xy_y_z = buffer_1010_pdpp[986];

    auto g_z_0_x_0_x_xy_z_x = buffer_1010_pdpp[987];

    auto g_z_0_x_0_x_xy_z_y = buffer_1010_pdpp[988];

    auto g_z_0_x_0_x_xy_z_z = buffer_1010_pdpp[989];

    auto g_z_0_x_0_x_xz_x_x = buffer_1010_pdpp[990];

    auto g_z_0_x_0_x_xz_x_y = buffer_1010_pdpp[991];

    auto g_z_0_x_0_x_xz_x_z = buffer_1010_pdpp[992];

    auto g_z_0_x_0_x_xz_y_x = buffer_1010_pdpp[993];

    auto g_z_0_x_0_x_xz_y_y = buffer_1010_pdpp[994];

    auto g_z_0_x_0_x_xz_y_z = buffer_1010_pdpp[995];

    auto g_z_0_x_0_x_xz_z_x = buffer_1010_pdpp[996];

    auto g_z_0_x_0_x_xz_z_y = buffer_1010_pdpp[997];

    auto g_z_0_x_0_x_xz_z_z = buffer_1010_pdpp[998];

    auto g_z_0_x_0_x_yy_x_x = buffer_1010_pdpp[999];

    auto g_z_0_x_0_x_yy_x_y = buffer_1010_pdpp[1000];

    auto g_z_0_x_0_x_yy_x_z = buffer_1010_pdpp[1001];

    auto g_z_0_x_0_x_yy_y_x = buffer_1010_pdpp[1002];

    auto g_z_0_x_0_x_yy_y_y = buffer_1010_pdpp[1003];

    auto g_z_0_x_0_x_yy_y_z = buffer_1010_pdpp[1004];

    auto g_z_0_x_0_x_yy_z_x = buffer_1010_pdpp[1005];

    auto g_z_0_x_0_x_yy_z_y = buffer_1010_pdpp[1006];

    auto g_z_0_x_0_x_yy_z_z = buffer_1010_pdpp[1007];

    auto g_z_0_x_0_x_yz_x_x = buffer_1010_pdpp[1008];

    auto g_z_0_x_0_x_yz_x_y = buffer_1010_pdpp[1009];

    auto g_z_0_x_0_x_yz_x_z = buffer_1010_pdpp[1010];

    auto g_z_0_x_0_x_yz_y_x = buffer_1010_pdpp[1011];

    auto g_z_0_x_0_x_yz_y_y = buffer_1010_pdpp[1012];

    auto g_z_0_x_0_x_yz_y_z = buffer_1010_pdpp[1013];

    auto g_z_0_x_0_x_yz_z_x = buffer_1010_pdpp[1014];

    auto g_z_0_x_0_x_yz_z_y = buffer_1010_pdpp[1015];

    auto g_z_0_x_0_x_yz_z_z = buffer_1010_pdpp[1016];

    auto g_z_0_x_0_x_zz_x_x = buffer_1010_pdpp[1017];

    auto g_z_0_x_0_x_zz_x_y = buffer_1010_pdpp[1018];

    auto g_z_0_x_0_x_zz_x_z = buffer_1010_pdpp[1019];

    auto g_z_0_x_0_x_zz_y_x = buffer_1010_pdpp[1020];

    auto g_z_0_x_0_x_zz_y_y = buffer_1010_pdpp[1021];

    auto g_z_0_x_0_x_zz_y_z = buffer_1010_pdpp[1022];

    auto g_z_0_x_0_x_zz_z_x = buffer_1010_pdpp[1023];

    auto g_z_0_x_0_x_zz_z_y = buffer_1010_pdpp[1024];

    auto g_z_0_x_0_x_zz_z_z = buffer_1010_pdpp[1025];

    auto g_z_0_x_0_y_xx_x_x = buffer_1010_pdpp[1026];

    auto g_z_0_x_0_y_xx_x_y = buffer_1010_pdpp[1027];

    auto g_z_0_x_0_y_xx_x_z = buffer_1010_pdpp[1028];

    auto g_z_0_x_0_y_xx_y_x = buffer_1010_pdpp[1029];

    auto g_z_0_x_0_y_xx_y_y = buffer_1010_pdpp[1030];

    auto g_z_0_x_0_y_xx_y_z = buffer_1010_pdpp[1031];

    auto g_z_0_x_0_y_xx_z_x = buffer_1010_pdpp[1032];

    auto g_z_0_x_0_y_xx_z_y = buffer_1010_pdpp[1033];

    auto g_z_0_x_0_y_xx_z_z = buffer_1010_pdpp[1034];

    auto g_z_0_x_0_y_xy_x_x = buffer_1010_pdpp[1035];

    auto g_z_0_x_0_y_xy_x_y = buffer_1010_pdpp[1036];

    auto g_z_0_x_0_y_xy_x_z = buffer_1010_pdpp[1037];

    auto g_z_0_x_0_y_xy_y_x = buffer_1010_pdpp[1038];

    auto g_z_0_x_0_y_xy_y_y = buffer_1010_pdpp[1039];

    auto g_z_0_x_0_y_xy_y_z = buffer_1010_pdpp[1040];

    auto g_z_0_x_0_y_xy_z_x = buffer_1010_pdpp[1041];

    auto g_z_0_x_0_y_xy_z_y = buffer_1010_pdpp[1042];

    auto g_z_0_x_0_y_xy_z_z = buffer_1010_pdpp[1043];

    auto g_z_0_x_0_y_xz_x_x = buffer_1010_pdpp[1044];

    auto g_z_0_x_0_y_xz_x_y = buffer_1010_pdpp[1045];

    auto g_z_0_x_0_y_xz_x_z = buffer_1010_pdpp[1046];

    auto g_z_0_x_0_y_xz_y_x = buffer_1010_pdpp[1047];

    auto g_z_0_x_0_y_xz_y_y = buffer_1010_pdpp[1048];

    auto g_z_0_x_0_y_xz_y_z = buffer_1010_pdpp[1049];

    auto g_z_0_x_0_y_xz_z_x = buffer_1010_pdpp[1050];

    auto g_z_0_x_0_y_xz_z_y = buffer_1010_pdpp[1051];

    auto g_z_0_x_0_y_xz_z_z = buffer_1010_pdpp[1052];

    auto g_z_0_x_0_y_yy_x_x = buffer_1010_pdpp[1053];

    auto g_z_0_x_0_y_yy_x_y = buffer_1010_pdpp[1054];

    auto g_z_0_x_0_y_yy_x_z = buffer_1010_pdpp[1055];

    auto g_z_0_x_0_y_yy_y_x = buffer_1010_pdpp[1056];

    auto g_z_0_x_0_y_yy_y_y = buffer_1010_pdpp[1057];

    auto g_z_0_x_0_y_yy_y_z = buffer_1010_pdpp[1058];

    auto g_z_0_x_0_y_yy_z_x = buffer_1010_pdpp[1059];

    auto g_z_0_x_0_y_yy_z_y = buffer_1010_pdpp[1060];

    auto g_z_0_x_0_y_yy_z_z = buffer_1010_pdpp[1061];

    auto g_z_0_x_0_y_yz_x_x = buffer_1010_pdpp[1062];

    auto g_z_0_x_0_y_yz_x_y = buffer_1010_pdpp[1063];

    auto g_z_0_x_0_y_yz_x_z = buffer_1010_pdpp[1064];

    auto g_z_0_x_0_y_yz_y_x = buffer_1010_pdpp[1065];

    auto g_z_0_x_0_y_yz_y_y = buffer_1010_pdpp[1066];

    auto g_z_0_x_0_y_yz_y_z = buffer_1010_pdpp[1067];

    auto g_z_0_x_0_y_yz_z_x = buffer_1010_pdpp[1068];

    auto g_z_0_x_0_y_yz_z_y = buffer_1010_pdpp[1069];

    auto g_z_0_x_0_y_yz_z_z = buffer_1010_pdpp[1070];

    auto g_z_0_x_0_y_zz_x_x = buffer_1010_pdpp[1071];

    auto g_z_0_x_0_y_zz_x_y = buffer_1010_pdpp[1072];

    auto g_z_0_x_0_y_zz_x_z = buffer_1010_pdpp[1073];

    auto g_z_0_x_0_y_zz_y_x = buffer_1010_pdpp[1074];

    auto g_z_0_x_0_y_zz_y_y = buffer_1010_pdpp[1075];

    auto g_z_0_x_0_y_zz_y_z = buffer_1010_pdpp[1076];

    auto g_z_0_x_0_y_zz_z_x = buffer_1010_pdpp[1077];

    auto g_z_0_x_0_y_zz_z_y = buffer_1010_pdpp[1078];

    auto g_z_0_x_0_y_zz_z_z = buffer_1010_pdpp[1079];

    auto g_z_0_x_0_z_xx_x_x = buffer_1010_pdpp[1080];

    auto g_z_0_x_0_z_xx_x_y = buffer_1010_pdpp[1081];

    auto g_z_0_x_0_z_xx_x_z = buffer_1010_pdpp[1082];

    auto g_z_0_x_0_z_xx_y_x = buffer_1010_pdpp[1083];

    auto g_z_0_x_0_z_xx_y_y = buffer_1010_pdpp[1084];

    auto g_z_0_x_0_z_xx_y_z = buffer_1010_pdpp[1085];

    auto g_z_0_x_0_z_xx_z_x = buffer_1010_pdpp[1086];

    auto g_z_0_x_0_z_xx_z_y = buffer_1010_pdpp[1087];

    auto g_z_0_x_0_z_xx_z_z = buffer_1010_pdpp[1088];

    auto g_z_0_x_0_z_xy_x_x = buffer_1010_pdpp[1089];

    auto g_z_0_x_0_z_xy_x_y = buffer_1010_pdpp[1090];

    auto g_z_0_x_0_z_xy_x_z = buffer_1010_pdpp[1091];

    auto g_z_0_x_0_z_xy_y_x = buffer_1010_pdpp[1092];

    auto g_z_0_x_0_z_xy_y_y = buffer_1010_pdpp[1093];

    auto g_z_0_x_0_z_xy_y_z = buffer_1010_pdpp[1094];

    auto g_z_0_x_0_z_xy_z_x = buffer_1010_pdpp[1095];

    auto g_z_0_x_0_z_xy_z_y = buffer_1010_pdpp[1096];

    auto g_z_0_x_0_z_xy_z_z = buffer_1010_pdpp[1097];

    auto g_z_0_x_0_z_xz_x_x = buffer_1010_pdpp[1098];

    auto g_z_0_x_0_z_xz_x_y = buffer_1010_pdpp[1099];

    auto g_z_0_x_0_z_xz_x_z = buffer_1010_pdpp[1100];

    auto g_z_0_x_0_z_xz_y_x = buffer_1010_pdpp[1101];

    auto g_z_0_x_0_z_xz_y_y = buffer_1010_pdpp[1102];

    auto g_z_0_x_0_z_xz_y_z = buffer_1010_pdpp[1103];

    auto g_z_0_x_0_z_xz_z_x = buffer_1010_pdpp[1104];

    auto g_z_0_x_0_z_xz_z_y = buffer_1010_pdpp[1105];

    auto g_z_0_x_0_z_xz_z_z = buffer_1010_pdpp[1106];

    auto g_z_0_x_0_z_yy_x_x = buffer_1010_pdpp[1107];

    auto g_z_0_x_0_z_yy_x_y = buffer_1010_pdpp[1108];

    auto g_z_0_x_0_z_yy_x_z = buffer_1010_pdpp[1109];

    auto g_z_0_x_0_z_yy_y_x = buffer_1010_pdpp[1110];

    auto g_z_0_x_0_z_yy_y_y = buffer_1010_pdpp[1111];

    auto g_z_0_x_0_z_yy_y_z = buffer_1010_pdpp[1112];

    auto g_z_0_x_0_z_yy_z_x = buffer_1010_pdpp[1113];

    auto g_z_0_x_0_z_yy_z_y = buffer_1010_pdpp[1114];

    auto g_z_0_x_0_z_yy_z_z = buffer_1010_pdpp[1115];

    auto g_z_0_x_0_z_yz_x_x = buffer_1010_pdpp[1116];

    auto g_z_0_x_0_z_yz_x_y = buffer_1010_pdpp[1117];

    auto g_z_0_x_0_z_yz_x_z = buffer_1010_pdpp[1118];

    auto g_z_0_x_0_z_yz_y_x = buffer_1010_pdpp[1119];

    auto g_z_0_x_0_z_yz_y_y = buffer_1010_pdpp[1120];

    auto g_z_0_x_0_z_yz_y_z = buffer_1010_pdpp[1121];

    auto g_z_0_x_0_z_yz_z_x = buffer_1010_pdpp[1122];

    auto g_z_0_x_0_z_yz_z_y = buffer_1010_pdpp[1123];

    auto g_z_0_x_0_z_yz_z_z = buffer_1010_pdpp[1124];

    auto g_z_0_x_0_z_zz_x_x = buffer_1010_pdpp[1125];

    auto g_z_0_x_0_z_zz_x_y = buffer_1010_pdpp[1126];

    auto g_z_0_x_0_z_zz_x_z = buffer_1010_pdpp[1127];

    auto g_z_0_x_0_z_zz_y_x = buffer_1010_pdpp[1128];

    auto g_z_0_x_0_z_zz_y_y = buffer_1010_pdpp[1129];

    auto g_z_0_x_0_z_zz_y_z = buffer_1010_pdpp[1130];

    auto g_z_0_x_0_z_zz_z_x = buffer_1010_pdpp[1131];

    auto g_z_0_x_0_z_zz_z_y = buffer_1010_pdpp[1132];

    auto g_z_0_x_0_z_zz_z_z = buffer_1010_pdpp[1133];

    auto g_z_0_y_0_x_xx_x_x = buffer_1010_pdpp[1134];

    auto g_z_0_y_0_x_xx_x_y = buffer_1010_pdpp[1135];

    auto g_z_0_y_0_x_xx_x_z = buffer_1010_pdpp[1136];

    auto g_z_0_y_0_x_xx_y_x = buffer_1010_pdpp[1137];

    auto g_z_0_y_0_x_xx_y_y = buffer_1010_pdpp[1138];

    auto g_z_0_y_0_x_xx_y_z = buffer_1010_pdpp[1139];

    auto g_z_0_y_0_x_xx_z_x = buffer_1010_pdpp[1140];

    auto g_z_0_y_0_x_xx_z_y = buffer_1010_pdpp[1141];

    auto g_z_0_y_0_x_xx_z_z = buffer_1010_pdpp[1142];

    auto g_z_0_y_0_x_xy_x_x = buffer_1010_pdpp[1143];

    auto g_z_0_y_0_x_xy_x_y = buffer_1010_pdpp[1144];

    auto g_z_0_y_0_x_xy_x_z = buffer_1010_pdpp[1145];

    auto g_z_0_y_0_x_xy_y_x = buffer_1010_pdpp[1146];

    auto g_z_0_y_0_x_xy_y_y = buffer_1010_pdpp[1147];

    auto g_z_0_y_0_x_xy_y_z = buffer_1010_pdpp[1148];

    auto g_z_0_y_0_x_xy_z_x = buffer_1010_pdpp[1149];

    auto g_z_0_y_0_x_xy_z_y = buffer_1010_pdpp[1150];

    auto g_z_0_y_0_x_xy_z_z = buffer_1010_pdpp[1151];

    auto g_z_0_y_0_x_xz_x_x = buffer_1010_pdpp[1152];

    auto g_z_0_y_0_x_xz_x_y = buffer_1010_pdpp[1153];

    auto g_z_0_y_0_x_xz_x_z = buffer_1010_pdpp[1154];

    auto g_z_0_y_0_x_xz_y_x = buffer_1010_pdpp[1155];

    auto g_z_0_y_0_x_xz_y_y = buffer_1010_pdpp[1156];

    auto g_z_0_y_0_x_xz_y_z = buffer_1010_pdpp[1157];

    auto g_z_0_y_0_x_xz_z_x = buffer_1010_pdpp[1158];

    auto g_z_0_y_0_x_xz_z_y = buffer_1010_pdpp[1159];

    auto g_z_0_y_0_x_xz_z_z = buffer_1010_pdpp[1160];

    auto g_z_0_y_0_x_yy_x_x = buffer_1010_pdpp[1161];

    auto g_z_0_y_0_x_yy_x_y = buffer_1010_pdpp[1162];

    auto g_z_0_y_0_x_yy_x_z = buffer_1010_pdpp[1163];

    auto g_z_0_y_0_x_yy_y_x = buffer_1010_pdpp[1164];

    auto g_z_0_y_0_x_yy_y_y = buffer_1010_pdpp[1165];

    auto g_z_0_y_0_x_yy_y_z = buffer_1010_pdpp[1166];

    auto g_z_0_y_0_x_yy_z_x = buffer_1010_pdpp[1167];

    auto g_z_0_y_0_x_yy_z_y = buffer_1010_pdpp[1168];

    auto g_z_0_y_0_x_yy_z_z = buffer_1010_pdpp[1169];

    auto g_z_0_y_0_x_yz_x_x = buffer_1010_pdpp[1170];

    auto g_z_0_y_0_x_yz_x_y = buffer_1010_pdpp[1171];

    auto g_z_0_y_0_x_yz_x_z = buffer_1010_pdpp[1172];

    auto g_z_0_y_0_x_yz_y_x = buffer_1010_pdpp[1173];

    auto g_z_0_y_0_x_yz_y_y = buffer_1010_pdpp[1174];

    auto g_z_0_y_0_x_yz_y_z = buffer_1010_pdpp[1175];

    auto g_z_0_y_0_x_yz_z_x = buffer_1010_pdpp[1176];

    auto g_z_0_y_0_x_yz_z_y = buffer_1010_pdpp[1177];

    auto g_z_0_y_0_x_yz_z_z = buffer_1010_pdpp[1178];

    auto g_z_0_y_0_x_zz_x_x = buffer_1010_pdpp[1179];

    auto g_z_0_y_0_x_zz_x_y = buffer_1010_pdpp[1180];

    auto g_z_0_y_0_x_zz_x_z = buffer_1010_pdpp[1181];

    auto g_z_0_y_0_x_zz_y_x = buffer_1010_pdpp[1182];

    auto g_z_0_y_0_x_zz_y_y = buffer_1010_pdpp[1183];

    auto g_z_0_y_0_x_zz_y_z = buffer_1010_pdpp[1184];

    auto g_z_0_y_0_x_zz_z_x = buffer_1010_pdpp[1185];

    auto g_z_0_y_0_x_zz_z_y = buffer_1010_pdpp[1186];

    auto g_z_0_y_0_x_zz_z_z = buffer_1010_pdpp[1187];

    auto g_z_0_y_0_y_xx_x_x = buffer_1010_pdpp[1188];

    auto g_z_0_y_0_y_xx_x_y = buffer_1010_pdpp[1189];

    auto g_z_0_y_0_y_xx_x_z = buffer_1010_pdpp[1190];

    auto g_z_0_y_0_y_xx_y_x = buffer_1010_pdpp[1191];

    auto g_z_0_y_0_y_xx_y_y = buffer_1010_pdpp[1192];

    auto g_z_0_y_0_y_xx_y_z = buffer_1010_pdpp[1193];

    auto g_z_0_y_0_y_xx_z_x = buffer_1010_pdpp[1194];

    auto g_z_0_y_0_y_xx_z_y = buffer_1010_pdpp[1195];

    auto g_z_0_y_0_y_xx_z_z = buffer_1010_pdpp[1196];

    auto g_z_0_y_0_y_xy_x_x = buffer_1010_pdpp[1197];

    auto g_z_0_y_0_y_xy_x_y = buffer_1010_pdpp[1198];

    auto g_z_0_y_0_y_xy_x_z = buffer_1010_pdpp[1199];

    auto g_z_0_y_0_y_xy_y_x = buffer_1010_pdpp[1200];

    auto g_z_0_y_0_y_xy_y_y = buffer_1010_pdpp[1201];

    auto g_z_0_y_0_y_xy_y_z = buffer_1010_pdpp[1202];

    auto g_z_0_y_0_y_xy_z_x = buffer_1010_pdpp[1203];

    auto g_z_0_y_0_y_xy_z_y = buffer_1010_pdpp[1204];

    auto g_z_0_y_0_y_xy_z_z = buffer_1010_pdpp[1205];

    auto g_z_0_y_0_y_xz_x_x = buffer_1010_pdpp[1206];

    auto g_z_0_y_0_y_xz_x_y = buffer_1010_pdpp[1207];

    auto g_z_0_y_0_y_xz_x_z = buffer_1010_pdpp[1208];

    auto g_z_0_y_0_y_xz_y_x = buffer_1010_pdpp[1209];

    auto g_z_0_y_0_y_xz_y_y = buffer_1010_pdpp[1210];

    auto g_z_0_y_0_y_xz_y_z = buffer_1010_pdpp[1211];

    auto g_z_0_y_0_y_xz_z_x = buffer_1010_pdpp[1212];

    auto g_z_0_y_0_y_xz_z_y = buffer_1010_pdpp[1213];

    auto g_z_0_y_0_y_xz_z_z = buffer_1010_pdpp[1214];

    auto g_z_0_y_0_y_yy_x_x = buffer_1010_pdpp[1215];

    auto g_z_0_y_0_y_yy_x_y = buffer_1010_pdpp[1216];

    auto g_z_0_y_0_y_yy_x_z = buffer_1010_pdpp[1217];

    auto g_z_0_y_0_y_yy_y_x = buffer_1010_pdpp[1218];

    auto g_z_0_y_0_y_yy_y_y = buffer_1010_pdpp[1219];

    auto g_z_0_y_0_y_yy_y_z = buffer_1010_pdpp[1220];

    auto g_z_0_y_0_y_yy_z_x = buffer_1010_pdpp[1221];

    auto g_z_0_y_0_y_yy_z_y = buffer_1010_pdpp[1222];

    auto g_z_0_y_0_y_yy_z_z = buffer_1010_pdpp[1223];

    auto g_z_0_y_0_y_yz_x_x = buffer_1010_pdpp[1224];

    auto g_z_0_y_0_y_yz_x_y = buffer_1010_pdpp[1225];

    auto g_z_0_y_0_y_yz_x_z = buffer_1010_pdpp[1226];

    auto g_z_0_y_0_y_yz_y_x = buffer_1010_pdpp[1227];

    auto g_z_0_y_0_y_yz_y_y = buffer_1010_pdpp[1228];

    auto g_z_0_y_0_y_yz_y_z = buffer_1010_pdpp[1229];

    auto g_z_0_y_0_y_yz_z_x = buffer_1010_pdpp[1230];

    auto g_z_0_y_0_y_yz_z_y = buffer_1010_pdpp[1231];

    auto g_z_0_y_0_y_yz_z_z = buffer_1010_pdpp[1232];

    auto g_z_0_y_0_y_zz_x_x = buffer_1010_pdpp[1233];

    auto g_z_0_y_0_y_zz_x_y = buffer_1010_pdpp[1234];

    auto g_z_0_y_0_y_zz_x_z = buffer_1010_pdpp[1235];

    auto g_z_0_y_0_y_zz_y_x = buffer_1010_pdpp[1236];

    auto g_z_0_y_0_y_zz_y_y = buffer_1010_pdpp[1237];

    auto g_z_0_y_0_y_zz_y_z = buffer_1010_pdpp[1238];

    auto g_z_0_y_0_y_zz_z_x = buffer_1010_pdpp[1239];

    auto g_z_0_y_0_y_zz_z_y = buffer_1010_pdpp[1240];

    auto g_z_0_y_0_y_zz_z_z = buffer_1010_pdpp[1241];

    auto g_z_0_y_0_z_xx_x_x = buffer_1010_pdpp[1242];

    auto g_z_0_y_0_z_xx_x_y = buffer_1010_pdpp[1243];

    auto g_z_0_y_0_z_xx_x_z = buffer_1010_pdpp[1244];

    auto g_z_0_y_0_z_xx_y_x = buffer_1010_pdpp[1245];

    auto g_z_0_y_0_z_xx_y_y = buffer_1010_pdpp[1246];

    auto g_z_0_y_0_z_xx_y_z = buffer_1010_pdpp[1247];

    auto g_z_0_y_0_z_xx_z_x = buffer_1010_pdpp[1248];

    auto g_z_0_y_0_z_xx_z_y = buffer_1010_pdpp[1249];

    auto g_z_0_y_0_z_xx_z_z = buffer_1010_pdpp[1250];

    auto g_z_0_y_0_z_xy_x_x = buffer_1010_pdpp[1251];

    auto g_z_0_y_0_z_xy_x_y = buffer_1010_pdpp[1252];

    auto g_z_0_y_0_z_xy_x_z = buffer_1010_pdpp[1253];

    auto g_z_0_y_0_z_xy_y_x = buffer_1010_pdpp[1254];

    auto g_z_0_y_0_z_xy_y_y = buffer_1010_pdpp[1255];

    auto g_z_0_y_0_z_xy_y_z = buffer_1010_pdpp[1256];

    auto g_z_0_y_0_z_xy_z_x = buffer_1010_pdpp[1257];

    auto g_z_0_y_0_z_xy_z_y = buffer_1010_pdpp[1258];

    auto g_z_0_y_0_z_xy_z_z = buffer_1010_pdpp[1259];

    auto g_z_0_y_0_z_xz_x_x = buffer_1010_pdpp[1260];

    auto g_z_0_y_0_z_xz_x_y = buffer_1010_pdpp[1261];

    auto g_z_0_y_0_z_xz_x_z = buffer_1010_pdpp[1262];

    auto g_z_0_y_0_z_xz_y_x = buffer_1010_pdpp[1263];

    auto g_z_0_y_0_z_xz_y_y = buffer_1010_pdpp[1264];

    auto g_z_0_y_0_z_xz_y_z = buffer_1010_pdpp[1265];

    auto g_z_0_y_0_z_xz_z_x = buffer_1010_pdpp[1266];

    auto g_z_0_y_0_z_xz_z_y = buffer_1010_pdpp[1267];

    auto g_z_0_y_0_z_xz_z_z = buffer_1010_pdpp[1268];

    auto g_z_0_y_0_z_yy_x_x = buffer_1010_pdpp[1269];

    auto g_z_0_y_0_z_yy_x_y = buffer_1010_pdpp[1270];

    auto g_z_0_y_0_z_yy_x_z = buffer_1010_pdpp[1271];

    auto g_z_0_y_0_z_yy_y_x = buffer_1010_pdpp[1272];

    auto g_z_0_y_0_z_yy_y_y = buffer_1010_pdpp[1273];

    auto g_z_0_y_0_z_yy_y_z = buffer_1010_pdpp[1274];

    auto g_z_0_y_0_z_yy_z_x = buffer_1010_pdpp[1275];

    auto g_z_0_y_0_z_yy_z_y = buffer_1010_pdpp[1276];

    auto g_z_0_y_0_z_yy_z_z = buffer_1010_pdpp[1277];

    auto g_z_0_y_0_z_yz_x_x = buffer_1010_pdpp[1278];

    auto g_z_0_y_0_z_yz_x_y = buffer_1010_pdpp[1279];

    auto g_z_0_y_0_z_yz_x_z = buffer_1010_pdpp[1280];

    auto g_z_0_y_0_z_yz_y_x = buffer_1010_pdpp[1281];

    auto g_z_0_y_0_z_yz_y_y = buffer_1010_pdpp[1282];

    auto g_z_0_y_0_z_yz_y_z = buffer_1010_pdpp[1283];

    auto g_z_0_y_0_z_yz_z_x = buffer_1010_pdpp[1284];

    auto g_z_0_y_0_z_yz_z_y = buffer_1010_pdpp[1285];

    auto g_z_0_y_0_z_yz_z_z = buffer_1010_pdpp[1286];

    auto g_z_0_y_0_z_zz_x_x = buffer_1010_pdpp[1287];

    auto g_z_0_y_0_z_zz_x_y = buffer_1010_pdpp[1288];

    auto g_z_0_y_0_z_zz_x_z = buffer_1010_pdpp[1289];

    auto g_z_0_y_0_z_zz_y_x = buffer_1010_pdpp[1290];

    auto g_z_0_y_0_z_zz_y_y = buffer_1010_pdpp[1291];

    auto g_z_0_y_0_z_zz_y_z = buffer_1010_pdpp[1292];

    auto g_z_0_y_0_z_zz_z_x = buffer_1010_pdpp[1293];

    auto g_z_0_y_0_z_zz_z_y = buffer_1010_pdpp[1294];

    auto g_z_0_y_0_z_zz_z_z = buffer_1010_pdpp[1295];

    auto g_z_0_z_0_x_xx_x_x = buffer_1010_pdpp[1296];

    auto g_z_0_z_0_x_xx_x_y = buffer_1010_pdpp[1297];

    auto g_z_0_z_0_x_xx_x_z = buffer_1010_pdpp[1298];

    auto g_z_0_z_0_x_xx_y_x = buffer_1010_pdpp[1299];

    auto g_z_0_z_0_x_xx_y_y = buffer_1010_pdpp[1300];

    auto g_z_0_z_0_x_xx_y_z = buffer_1010_pdpp[1301];

    auto g_z_0_z_0_x_xx_z_x = buffer_1010_pdpp[1302];

    auto g_z_0_z_0_x_xx_z_y = buffer_1010_pdpp[1303];

    auto g_z_0_z_0_x_xx_z_z = buffer_1010_pdpp[1304];

    auto g_z_0_z_0_x_xy_x_x = buffer_1010_pdpp[1305];

    auto g_z_0_z_0_x_xy_x_y = buffer_1010_pdpp[1306];

    auto g_z_0_z_0_x_xy_x_z = buffer_1010_pdpp[1307];

    auto g_z_0_z_0_x_xy_y_x = buffer_1010_pdpp[1308];

    auto g_z_0_z_0_x_xy_y_y = buffer_1010_pdpp[1309];

    auto g_z_0_z_0_x_xy_y_z = buffer_1010_pdpp[1310];

    auto g_z_0_z_0_x_xy_z_x = buffer_1010_pdpp[1311];

    auto g_z_0_z_0_x_xy_z_y = buffer_1010_pdpp[1312];

    auto g_z_0_z_0_x_xy_z_z = buffer_1010_pdpp[1313];

    auto g_z_0_z_0_x_xz_x_x = buffer_1010_pdpp[1314];

    auto g_z_0_z_0_x_xz_x_y = buffer_1010_pdpp[1315];

    auto g_z_0_z_0_x_xz_x_z = buffer_1010_pdpp[1316];

    auto g_z_0_z_0_x_xz_y_x = buffer_1010_pdpp[1317];

    auto g_z_0_z_0_x_xz_y_y = buffer_1010_pdpp[1318];

    auto g_z_0_z_0_x_xz_y_z = buffer_1010_pdpp[1319];

    auto g_z_0_z_0_x_xz_z_x = buffer_1010_pdpp[1320];

    auto g_z_0_z_0_x_xz_z_y = buffer_1010_pdpp[1321];

    auto g_z_0_z_0_x_xz_z_z = buffer_1010_pdpp[1322];

    auto g_z_0_z_0_x_yy_x_x = buffer_1010_pdpp[1323];

    auto g_z_0_z_0_x_yy_x_y = buffer_1010_pdpp[1324];

    auto g_z_0_z_0_x_yy_x_z = buffer_1010_pdpp[1325];

    auto g_z_0_z_0_x_yy_y_x = buffer_1010_pdpp[1326];

    auto g_z_0_z_0_x_yy_y_y = buffer_1010_pdpp[1327];

    auto g_z_0_z_0_x_yy_y_z = buffer_1010_pdpp[1328];

    auto g_z_0_z_0_x_yy_z_x = buffer_1010_pdpp[1329];

    auto g_z_0_z_0_x_yy_z_y = buffer_1010_pdpp[1330];

    auto g_z_0_z_0_x_yy_z_z = buffer_1010_pdpp[1331];

    auto g_z_0_z_0_x_yz_x_x = buffer_1010_pdpp[1332];

    auto g_z_0_z_0_x_yz_x_y = buffer_1010_pdpp[1333];

    auto g_z_0_z_0_x_yz_x_z = buffer_1010_pdpp[1334];

    auto g_z_0_z_0_x_yz_y_x = buffer_1010_pdpp[1335];

    auto g_z_0_z_0_x_yz_y_y = buffer_1010_pdpp[1336];

    auto g_z_0_z_0_x_yz_y_z = buffer_1010_pdpp[1337];

    auto g_z_0_z_0_x_yz_z_x = buffer_1010_pdpp[1338];

    auto g_z_0_z_0_x_yz_z_y = buffer_1010_pdpp[1339];

    auto g_z_0_z_0_x_yz_z_z = buffer_1010_pdpp[1340];

    auto g_z_0_z_0_x_zz_x_x = buffer_1010_pdpp[1341];

    auto g_z_0_z_0_x_zz_x_y = buffer_1010_pdpp[1342];

    auto g_z_0_z_0_x_zz_x_z = buffer_1010_pdpp[1343];

    auto g_z_0_z_0_x_zz_y_x = buffer_1010_pdpp[1344];

    auto g_z_0_z_0_x_zz_y_y = buffer_1010_pdpp[1345];

    auto g_z_0_z_0_x_zz_y_z = buffer_1010_pdpp[1346];

    auto g_z_0_z_0_x_zz_z_x = buffer_1010_pdpp[1347];

    auto g_z_0_z_0_x_zz_z_y = buffer_1010_pdpp[1348];

    auto g_z_0_z_0_x_zz_z_z = buffer_1010_pdpp[1349];

    auto g_z_0_z_0_y_xx_x_x = buffer_1010_pdpp[1350];

    auto g_z_0_z_0_y_xx_x_y = buffer_1010_pdpp[1351];

    auto g_z_0_z_0_y_xx_x_z = buffer_1010_pdpp[1352];

    auto g_z_0_z_0_y_xx_y_x = buffer_1010_pdpp[1353];

    auto g_z_0_z_0_y_xx_y_y = buffer_1010_pdpp[1354];

    auto g_z_0_z_0_y_xx_y_z = buffer_1010_pdpp[1355];

    auto g_z_0_z_0_y_xx_z_x = buffer_1010_pdpp[1356];

    auto g_z_0_z_0_y_xx_z_y = buffer_1010_pdpp[1357];

    auto g_z_0_z_0_y_xx_z_z = buffer_1010_pdpp[1358];

    auto g_z_0_z_0_y_xy_x_x = buffer_1010_pdpp[1359];

    auto g_z_0_z_0_y_xy_x_y = buffer_1010_pdpp[1360];

    auto g_z_0_z_0_y_xy_x_z = buffer_1010_pdpp[1361];

    auto g_z_0_z_0_y_xy_y_x = buffer_1010_pdpp[1362];

    auto g_z_0_z_0_y_xy_y_y = buffer_1010_pdpp[1363];

    auto g_z_0_z_0_y_xy_y_z = buffer_1010_pdpp[1364];

    auto g_z_0_z_0_y_xy_z_x = buffer_1010_pdpp[1365];

    auto g_z_0_z_0_y_xy_z_y = buffer_1010_pdpp[1366];

    auto g_z_0_z_0_y_xy_z_z = buffer_1010_pdpp[1367];

    auto g_z_0_z_0_y_xz_x_x = buffer_1010_pdpp[1368];

    auto g_z_0_z_0_y_xz_x_y = buffer_1010_pdpp[1369];

    auto g_z_0_z_0_y_xz_x_z = buffer_1010_pdpp[1370];

    auto g_z_0_z_0_y_xz_y_x = buffer_1010_pdpp[1371];

    auto g_z_0_z_0_y_xz_y_y = buffer_1010_pdpp[1372];

    auto g_z_0_z_0_y_xz_y_z = buffer_1010_pdpp[1373];

    auto g_z_0_z_0_y_xz_z_x = buffer_1010_pdpp[1374];

    auto g_z_0_z_0_y_xz_z_y = buffer_1010_pdpp[1375];

    auto g_z_0_z_0_y_xz_z_z = buffer_1010_pdpp[1376];

    auto g_z_0_z_0_y_yy_x_x = buffer_1010_pdpp[1377];

    auto g_z_0_z_0_y_yy_x_y = buffer_1010_pdpp[1378];

    auto g_z_0_z_0_y_yy_x_z = buffer_1010_pdpp[1379];

    auto g_z_0_z_0_y_yy_y_x = buffer_1010_pdpp[1380];

    auto g_z_0_z_0_y_yy_y_y = buffer_1010_pdpp[1381];

    auto g_z_0_z_0_y_yy_y_z = buffer_1010_pdpp[1382];

    auto g_z_0_z_0_y_yy_z_x = buffer_1010_pdpp[1383];

    auto g_z_0_z_0_y_yy_z_y = buffer_1010_pdpp[1384];

    auto g_z_0_z_0_y_yy_z_z = buffer_1010_pdpp[1385];

    auto g_z_0_z_0_y_yz_x_x = buffer_1010_pdpp[1386];

    auto g_z_0_z_0_y_yz_x_y = buffer_1010_pdpp[1387];

    auto g_z_0_z_0_y_yz_x_z = buffer_1010_pdpp[1388];

    auto g_z_0_z_0_y_yz_y_x = buffer_1010_pdpp[1389];

    auto g_z_0_z_0_y_yz_y_y = buffer_1010_pdpp[1390];

    auto g_z_0_z_0_y_yz_y_z = buffer_1010_pdpp[1391];

    auto g_z_0_z_0_y_yz_z_x = buffer_1010_pdpp[1392];

    auto g_z_0_z_0_y_yz_z_y = buffer_1010_pdpp[1393];

    auto g_z_0_z_0_y_yz_z_z = buffer_1010_pdpp[1394];

    auto g_z_0_z_0_y_zz_x_x = buffer_1010_pdpp[1395];

    auto g_z_0_z_0_y_zz_x_y = buffer_1010_pdpp[1396];

    auto g_z_0_z_0_y_zz_x_z = buffer_1010_pdpp[1397];

    auto g_z_0_z_0_y_zz_y_x = buffer_1010_pdpp[1398];

    auto g_z_0_z_0_y_zz_y_y = buffer_1010_pdpp[1399];

    auto g_z_0_z_0_y_zz_y_z = buffer_1010_pdpp[1400];

    auto g_z_0_z_0_y_zz_z_x = buffer_1010_pdpp[1401];

    auto g_z_0_z_0_y_zz_z_y = buffer_1010_pdpp[1402];

    auto g_z_0_z_0_y_zz_z_z = buffer_1010_pdpp[1403];

    auto g_z_0_z_0_z_xx_x_x = buffer_1010_pdpp[1404];

    auto g_z_0_z_0_z_xx_x_y = buffer_1010_pdpp[1405];

    auto g_z_0_z_0_z_xx_x_z = buffer_1010_pdpp[1406];

    auto g_z_0_z_0_z_xx_y_x = buffer_1010_pdpp[1407];

    auto g_z_0_z_0_z_xx_y_y = buffer_1010_pdpp[1408];

    auto g_z_0_z_0_z_xx_y_z = buffer_1010_pdpp[1409];

    auto g_z_0_z_0_z_xx_z_x = buffer_1010_pdpp[1410];

    auto g_z_0_z_0_z_xx_z_y = buffer_1010_pdpp[1411];

    auto g_z_0_z_0_z_xx_z_z = buffer_1010_pdpp[1412];

    auto g_z_0_z_0_z_xy_x_x = buffer_1010_pdpp[1413];

    auto g_z_0_z_0_z_xy_x_y = buffer_1010_pdpp[1414];

    auto g_z_0_z_0_z_xy_x_z = buffer_1010_pdpp[1415];

    auto g_z_0_z_0_z_xy_y_x = buffer_1010_pdpp[1416];

    auto g_z_0_z_0_z_xy_y_y = buffer_1010_pdpp[1417];

    auto g_z_0_z_0_z_xy_y_z = buffer_1010_pdpp[1418];

    auto g_z_0_z_0_z_xy_z_x = buffer_1010_pdpp[1419];

    auto g_z_0_z_0_z_xy_z_y = buffer_1010_pdpp[1420];

    auto g_z_0_z_0_z_xy_z_z = buffer_1010_pdpp[1421];

    auto g_z_0_z_0_z_xz_x_x = buffer_1010_pdpp[1422];

    auto g_z_0_z_0_z_xz_x_y = buffer_1010_pdpp[1423];

    auto g_z_0_z_0_z_xz_x_z = buffer_1010_pdpp[1424];

    auto g_z_0_z_0_z_xz_y_x = buffer_1010_pdpp[1425];

    auto g_z_0_z_0_z_xz_y_y = buffer_1010_pdpp[1426];

    auto g_z_0_z_0_z_xz_y_z = buffer_1010_pdpp[1427];

    auto g_z_0_z_0_z_xz_z_x = buffer_1010_pdpp[1428];

    auto g_z_0_z_0_z_xz_z_y = buffer_1010_pdpp[1429];

    auto g_z_0_z_0_z_xz_z_z = buffer_1010_pdpp[1430];

    auto g_z_0_z_0_z_yy_x_x = buffer_1010_pdpp[1431];

    auto g_z_0_z_0_z_yy_x_y = buffer_1010_pdpp[1432];

    auto g_z_0_z_0_z_yy_x_z = buffer_1010_pdpp[1433];

    auto g_z_0_z_0_z_yy_y_x = buffer_1010_pdpp[1434];

    auto g_z_0_z_0_z_yy_y_y = buffer_1010_pdpp[1435];

    auto g_z_0_z_0_z_yy_y_z = buffer_1010_pdpp[1436];

    auto g_z_0_z_0_z_yy_z_x = buffer_1010_pdpp[1437];

    auto g_z_0_z_0_z_yy_z_y = buffer_1010_pdpp[1438];

    auto g_z_0_z_0_z_yy_z_z = buffer_1010_pdpp[1439];

    auto g_z_0_z_0_z_yz_x_x = buffer_1010_pdpp[1440];

    auto g_z_0_z_0_z_yz_x_y = buffer_1010_pdpp[1441];

    auto g_z_0_z_0_z_yz_x_z = buffer_1010_pdpp[1442];

    auto g_z_0_z_0_z_yz_y_x = buffer_1010_pdpp[1443];

    auto g_z_0_z_0_z_yz_y_y = buffer_1010_pdpp[1444];

    auto g_z_0_z_0_z_yz_y_z = buffer_1010_pdpp[1445];

    auto g_z_0_z_0_z_yz_z_x = buffer_1010_pdpp[1446];

    auto g_z_0_z_0_z_yz_z_y = buffer_1010_pdpp[1447];

    auto g_z_0_z_0_z_yz_z_z = buffer_1010_pdpp[1448];

    auto g_z_0_z_0_z_zz_x_x = buffer_1010_pdpp[1449];

    auto g_z_0_z_0_z_zz_x_y = buffer_1010_pdpp[1450];

    auto g_z_0_z_0_z_zz_x_z = buffer_1010_pdpp[1451];

    auto g_z_0_z_0_z_zz_y_x = buffer_1010_pdpp[1452];

    auto g_z_0_z_0_z_zz_y_y = buffer_1010_pdpp[1453];

    auto g_z_0_z_0_z_zz_y_z = buffer_1010_pdpp[1454];

    auto g_z_0_z_0_z_zz_z_x = buffer_1010_pdpp[1455];

    auto g_z_0_z_0_z_zz_z_y = buffer_1010_pdpp[1456];

    auto g_z_0_z_0_z_zz_z_z = buffer_1010_pdpp[1457];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_0_xx_xx_x, g_0_xx_xx_y, g_0_xx_xx_z, g_x_0_x_0_x_xx_x_x, g_x_0_x_0_x_xx_x_y, g_x_0_x_0_x_xx_x_z, g_xx_xx_0_x, g_xx_xx_0_y, g_xx_xx_0_z, g_xx_xx_xx_x, g_xx_xx_xx_y, g_xx_xx_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xx_x_x[i] = g_0_xx_0_x[i] - 2.0 * g_0_xx_xx_x[i] * c_exps[i] - 2.0 * g_xx_xx_0_x[i] * a_exp + 4.0 * g_xx_xx_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_x_y[i] = g_0_xx_0_y[i] - 2.0 * g_0_xx_xx_y[i] * c_exps[i] - 2.0 * g_xx_xx_0_y[i] * a_exp + 4.0 * g_xx_xx_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_x_z[i] = g_0_xx_0_z[i] - 2.0 * g_0_xx_xx_z[i] * c_exps[i] - 2.0 * g_xx_xx_0_z[i] * a_exp + 4.0 * g_xx_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_xx_xy_x, g_0_xx_xy_y, g_0_xx_xy_z, g_x_0_x_0_x_xx_y_x, g_x_0_x_0_x_xx_y_y, g_x_0_x_0_x_xx_y_z, g_xx_xx_xy_x, g_xx_xx_xy_y, g_xx_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xx_y_x[i] = -2.0 * g_0_xx_xy_x[i] * c_exps[i] + 4.0 * g_xx_xx_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_y_y[i] = -2.0 * g_0_xx_xy_y[i] * c_exps[i] + 4.0 * g_xx_xx_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_y_z[i] = -2.0 * g_0_xx_xy_z[i] * c_exps[i] + 4.0 * g_xx_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_xx_xz_x, g_0_xx_xz_y, g_0_xx_xz_z, g_x_0_x_0_x_xx_z_x, g_x_0_x_0_x_xx_z_y, g_x_0_x_0_x_xx_z_z, g_xx_xx_xz_x, g_xx_xx_xz_y, g_xx_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xx_z_x[i] = -2.0 * g_0_xx_xz_x[i] * c_exps[i] + 4.0 * g_xx_xx_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_z_y[i] = -2.0 * g_0_xx_xz_y[i] * c_exps[i] + 4.0 * g_xx_xx_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_z_z[i] = -2.0 * g_0_xx_xz_z[i] * c_exps[i] + 4.0 * g_xx_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_0_xy_xx_x, g_0_xy_xx_y, g_0_xy_xx_z, g_x_0_x_0_x_xy_x_x, g_x_0_x_0_x_xy_x_y, g_x_0_x_0_x_xy_x_z, g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z, g_xx_xy_xx_x, g_xx_xy_xx_y, g_xx_xy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xy_x_x[i] = g_0_xy_0_x[i] - 2.0 * g_0_xy_xx_x[i] * c_exps[i] - 2.0 * g_xx_xy_0_x[i] * a_exp + 4.0 * g_xx_xy_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_x_y[i] = g_0_xy_0_y[i] - 2.0 * g_0_xy_xx_y[i] * c_exps[i] - 2.0 * g_xx_xy_0_y[i] * a_exp + 4.0 * g_xx_xy_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_x_z[i] = g_0_xy_0_z[i] - 2.0 * g_0_xy_xx_z[i] * c_exps[i] - 2.0 * g_xx_xy_0_z[i] * a_exp + 4.0 * g_xx_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_xy_xy_x, g_0_xy_xy_y, g_0_xy_xy_z, g_x_0_x_0_x_xy_y_x, g_x_0_x_0_x_xy_y_y, g_x_0_x_0_x_xy_y_z, g_xx_xy_xy_x, g_xx_xy_xy_y, g_xx_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xy_y_x[i] = -2.0 * g_0_xy_xy_x[i] * c_exps[i] + 4.0 * g_xx_xy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_y_y[i] = -2.0 * g_0_xy_xy_y[i] * c_exps[i] + 4.0 * g_xx_xy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_y_z[i] = -2.0 * g_0_xy_xy_z[i] * c_exps[i] + 4.0 * g_xx_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_xy_xz_x, g_0_xy_xz_y, g_0_xy_xz_z, g_x_0_x_0_x_xy_z_x, g_x_0_x_0_x_xy_z_y, g_x_0_x_0_x_xy_z_z, g_xx_xy_xz_x, g_xx_xy_xz_y, g_xx_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xy_z_x[i] = -2.0 * g_0_xy_xz_x[i] * c_exps[i] + 4.0 * g_xx_xy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_z_y[i] = -2.0 * g_0_xy_xz_y[i] * c_exps[i] + 4.0 * g_xx_xy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_z_z[i] = -2.0 * g_0_xy_xz_z[i] * c_exps[i] + 4.0 * g_xx_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_0_xz_xx_x, g_0_xz_xx_y, g_0_xz_xx_z, g_x_0_x_0_x_xz_x_x, g_x_0_x_0_x_xz_x_y, g_x_0_x_0_x_xz_x_z, g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z, g_xx_xz_xx_x, g_xx_xz_xx_y, g_xx_xz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xz_x_x[i] = g_0_xz_0_x[i] - 2.0 * g_0_xz_xx_x[i] * c_exps[i] - 2.0 * g_xx_xz_0_x[i] * a_exp + 4.0 * g_xx_xz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_x_y[i] = g_0_xz_0_y[i] - 2.0 * g_0_xz_xx_y[i] * c_exps[i] - 2.0 * g_xx_xz_0_y[i] * a_exp + 4.0 * g_xx_xz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_x_z[i] = g_0_xz_0_z[i] - 2.0 * g_0_xz_xx_z[i] * c_exps[i] - 2.0 * g_xx_xz_0_z[i] * a_exp + 4.0 * g_xx_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_0_xz_xy_x, g_0_xz_xy_y, g_0_xz_xy_z, g_x_0_x_0_x_xz_y_x, g_x_0_x_0_x_xz_y_y, g_x_0_x_0_x_xz_y_z, g_xx_xz_xy_x, g_xx_xz_xy_y, g_xx_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xz_y_x[i] = -2.0 * g_0_xz_xy_x[i] * c_exps[i] + 4.0 * g_xx_xz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_y_y[i] = -2.0 * g_0_xz_xy_y[i] * c_exps[i] + 4.0 * g_xx_xz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_y_z[i] = -2.0 * g_0_xz_xy_z[i] * c_exps[i] + 4.0 * g_xx_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_0_xz_xz_x, g_0_xz_xz_y, g_0_xz_xz_z, g_x_0_x_0_x_xz_z_x, g_x_0_x_0_x_xz_z_y, g_x_0_x_0_x_xz_z_z, g_xx_xz_xz_x, g_xx_xz_xz_y, g_xx_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xz_z_x[i] = -2.0 * g_0_xz_xz_x[i] * c_exps[i] + 4.0 * g_xx_xz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_z_y[i] = -2.0 * g_0_xz_xz_y[i] * c_exps[i] + 4.0 * g_xx_xz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_z_z[i] = -2.0 * g_0_xz_xz_z[i] * c_exps[i] + 4.0 * g_xx_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_0_yy_xx_x, g_0_yy_xx_y, g_0_yy_xx_z, g_x_0_x_0_x_yy_x_x, g_x_0_x_0_x_yy_x_y, g_x_0_x_0_x_yy_x_z, g_xx_yy_0_x, g_xx_yy_0_y, g_xx_yy_0_z, g_xx_yy_xx_x, g_xx_yy_xx_y, g_xx_yy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_yy_x_x[i] = g_0_yy_0_x[i] - 2.0 * g_0_yy_xx_x[i] * c_exps[i] - 2.0 * g_xx_yy_0_x[i] * a_exp + 4.0 * g_xx_yy_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_x_y[i] = g_0_yy_0_y[i] - 2.0 * g_0_yy_xx_y[i] * c_exps[i] - 2.0 * g_xx_yy_0_y[i] * a_exp + 4.0 * g_xx_yy_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_x_z[i] = g_0_yy_0_z[i] - 2.0 * g_0_yy_xx_z[i] * c_exps[i] - 2.0 * g_xx_yy_0_z[i] * a_exp + 4.0 * g_xx_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_0_yy_xy_x, g_0_yy_xy_y, g_0_yy_xy_z, g_x_0_x_0_x_yy_y_x, g_x_0_x_0_x_yy_y_y, g_x_0_x_0_x_yy_y_z, g_xx_yy_xy_x, g_xx_yy_xy_y, g_xx_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_yy_y_x[i] = -2.0 * g_0_yy_xy_x[i] * c_exps[i] + 4.0 * g_xx_yy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_y_y[i] = -2.0 * g_0_yy_xy_y[i] * c_exps[i] + 4.0 * g_xx_yy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_y_z[i] = -2.0 * g_0_yy_xy_z[i] * c_exps[i] + 4.0 * g_xx_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_0_yy_xz_x, g_0_yy_xz_y, g_0_yy_xz_z, g_x_0_x_0_x_yy_z_x, g_x_0_x_0_x_yy_z_y, g_x_0_x_0_x_yy_z_z, g_xx_yy_xz_x, g_xx_yy_xz_y, g_xx_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_yy_z_x[i] = -2.0 * g_0_yy_xz_x[i] * c_exps[i] + 4.0 * g_xx_yy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_z_y[i] = -2.0 * g_0_yy_xz_y[i] * c_exps[i] + 4.0 * g_xx_yy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_z_z[i] = -2.0 * g_0_yy_xz_z[i] * c_exps[i] + 4.0 * g_xx_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_0_yz_xx_x, g_0_yz_xx_y, g_0_yz_xx_z, g_x_0_x_0_x_yz_x_x, g_x_0_x_0_x_yz_x_y, g_x_0_x_0_x_yz_x_z, g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z, g_xx_yz_xx_x, g_xx_yz_xx_y, g_xx_yz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_yz_x_x[i] = g_0_yz_0_x[i] - 2.0 * g_0_yz_xx_x[i] * c_exps[i] - 2.0 * g_xx_yz_0_x[i] * a_exp + 4.0 * g_xx_yz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_x_y[i] = g_0_yz_0_y[i] - 2.0 * g_0_yz_xx_y[i] * c_exps[i] - 2.0 * g_xx_yz_0_y[i] * a_exp + 4.0 * g_xx_yz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_x_z[i] = g_0_yz_0_z[i] - 2.0 * g_0_yz_xx_z[i] * c_exps[i] - 2.0 * g_xx_yz_0_z[i] * a_exp + 4.0 * g_xx_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_0_yz_xy_x, g_0_yz_xy_y, g_0_yz_xy_z, g_x_0_x_0_x_yz_y_x, g_x_0_x_0_x_yz_y_y, g_x_0_x_0_x_yz_y_z, g_xx_yz_xy_x, g_xx_yz_xy_y, g_xx_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_yz_y_x[i] = -2.0 * g_0_yz_xy_x[i] * c_exps[i] + 4.0 * g_xx_yz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_y_y[i] = -2.0 * g_0_yz_xy_y[i] * c_exps[i] + 4.0 * g_xx_yz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_y_z[i] = -2.0 * g_0_yz_xy_z[i] * c_exps[i] + 4.0 * g_xx_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_0_yz_xz_x, g_0_yz_xz_y, g_0_yz_xz_z, g_x_0_x_0_x_yz_z_x, g_x_0_x_0_x_yz_z_y, g_x_0_x_0_x_yz_z_z, g_xx_yz_xz_x, g_xx_yz_xz_y, g_xx_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_yz_z_x[i] = -2.0 * g_0_yz_xz_x[i] * c_exps[i] + 4.0 * g_xx_yz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_z_y[i] = -2.0 * g_0_yz_xz_y[i] * c_exps[i] + 4.0 * g_xx_yz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_z_z[i] = -2.0 * g_0_yz_xz_z[i] * c_exps[i] + 4.0 * g_xx_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_0_zz_xx_x, g_0_zz_xx_y, g_0_zz_xx_z, g_x_0_x_0_x_zz_x_x, g_x_0_x_0_x_zz_x_y, g_x_0_x_0_x_zz_x_z, g_xx_zz_0_x, g_xx_zz_0_y, g_xx_zz_0_z, g_xx_zz_xx_x, g_xx_zz_xx_y, g_xx_zz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_zz_x_x[i] = g_0_zz_0_x[i] - 2.0 * g_0_zz_xx_x[i] * c_exps[i] - 2.0 * g_xx_zz_0_x[i] * a_exp + 4.0 * g_xx_zz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_x_y[i] = g_0_zz_0_y[i] - 2.0 * g_0_zz_xx_y[i] * c_exps[i] - 2.0 * g_xx_zz_0_y[i] * a_exp + 4.0 * g_xx_zz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_x_z[i] = g_0_zz_0_z[i] - 2.0 * g_0_zz_xx_z[i] * c_exps[i] - 2.0 * g_xx_zz_0_z[i] * a_exp + 4.0 * g_xx_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_0_zz_xy_x, g_0_zz_xy_y, g_0_zz_xy_z, g_x_0_x_0_x_zz_y_x, g_x_0_x_0_x_zz_y_y, g_x_0_x_0_x_zz_y_z, g_xx_zz_xy_x, g_xx_zz_xy_y, g_xx_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_zz_y_x[i] = -2.0 * g_0_zz_xy_x[i] * c_exps[i] + 4.0 * g_xx_zz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_y_y[i] = -2.0 * g_0_zz_xy_y[i] * c_exps[i] + 4.0 * g_xx_zz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_y_z[i] = -2.0 * g_0_zz_xy_z[i] * c_exps[i] + 4.0 * g_xx_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_0_zz_xz_x, g_0_zz_xz_y, g_0_zz_xz_z, g_x_0_x_0_x_zz_z_x, g_x_0_x_0_x_zz_z_y, g_x_0_x_0_x_zz_z_z, g_xx_zz_xz_x, g_xx_zz_xz_y, g_xx_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_zz_z_x[i] = -2.0 * g_0_zz_xz_x[i] * c_exps[i] + 4.0 * g_xx_zz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_z_y[i] = -2.0 * g_0_zz_xz_y[i] * c_exps[i] + 4.0 * g_xx_zz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_z_z[i] = -2.0 * g_0_zz_xz_z[i] * c_exps[i] + 4.0 * g_xx_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_0_x_0_y_xx_x_x, g_x_0_x_0_y_xx_x_y, g_x_0_x_0_y_xx_x_z, g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xy_xx_xx_x, g_xy_xx_xx_y, g_xy_xx_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xx_x_x[i] = -2.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xy_xx_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_x_y[i] = -2.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xy_xx_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_x_z[i] = -2.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xy_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_0_x_0_y_xx_y_x, g_x_0_x_0_y_xx_y_y, g_x_0_x_0_y_xx_y_z, g_xy_xx_xy_x, g_xy_xx_xy_y, g_xy_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xx_y_x[i] = 4.0 * g_xy_xx_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_y_y[i] = 4.0 * g_xy_xx_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_y_z[i] = 4.0 * g_xy_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_0_x_0_y_xx_z_x, g_x_0_x_0_y_xx_z_y, g_x_0_x_0_y_xx_z_z, g_xy_xx_xz_x, g_xy_xx_xz_y, g_xy_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xx_z_x[i] = 4.0 * g_xy_xx_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_z_y[i] = 4.0 * g_xy_xx_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_z_z[i] = 4.0 * g_xy_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_0_x_0_y_xy_x_x, g_x_0_x_0_y_xy_x_y, g_x_0_x_0_y_xy_x_z, g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xy_xy_xx_x, g_xy_xy_xx_y, g_xy_xy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xy_x_x[i] = -2.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xy_xy_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_x_y[i] = -2.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xy_xy_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_x_z[i] = -2.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xy_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_0_x_0_y_xy_y_x, g_x_0_x_0_y_xy_y_y, g_x_0_x_0_y_xy_y_z, g_xy_xy_xy_x, g_xy_xy_xy_y, g_xy_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xy_y_x[i] = 4.0 * g_xy_xy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_y_y[i] = 4.0 * g_xy_xy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_y_z[i] = 4.0 * g_xy_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_0_x_0_y_xy_z_x, g_x_0_x_0_y_xy_z_y, g_x_0_x_0_y_xy_z_z, g_xy_xy_xz_x, g_xy_xy_xz_y, g_xy_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xy_z_x[i] = 4.0 * g_xy_xy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_z_y[i] = 4.0 * g_xy_xy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_z_z[i] = 4.0 * g_xy_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_0_x_0_y_xz_x_x, g_x_0_x_0_y_xz_x_y, g_x_0_x_0_y_xz_x_z, g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xy_xz_xx_x, g_xy_xz_xx_y, g_xy_xz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xz_x_x[i] = -2.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xy_xz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_x_y[i] = -2.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xy_xz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_x_z[i] = -2.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xy_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_0_x_0_y_xz_y_x, g_x_0_x_0_y_xz_y_y, g_x_0_x_0_y_xz_y_z, g_xy_xz_xy_x, g_xy_xz_xy_y, g_xy_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xz_y_x[i] = 4.0 * g_xy_xz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_y_y[i] = 4.0 * g_xy_xz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_y_z[i] = 4.0 * g_xy_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_0_x_0_y_xz_z_x, g_x_0_x_0_y_xz_z_y, g_x_0_x_0_y_xz_z_z, g_xy_xz_xz_x, g_xy_xz_xz_y, g_xy_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xz_z_x[i] = 4.0 * g_xy_xz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_z_y[i] = 4.0 * g_xy_xz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_z_z[i] = 4.0 * g_xy_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_0_x_0_y_yy_x_x, g_x_0_x_0_y_yy_x_y, g_x_0_x_0_y_yy_x_z, g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xy_yy_xx_x, g_xy_yy_xx_y, g_xy_yy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_yy_x_x[i] = -2.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xy_yy_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_x_y[i] = -2.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xy_yy_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_x_z[i] = -2.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xy_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_0_x_0_y_yy_y_x, g_x_0_x_0_y_yy_y_y, g_x_0_x_0_y_yy_y_z, g_xy_yy_xy_x, g_xy_yy_xy_y, g_xy_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_yy_y_x[i] = 4.0 * g_xy_yy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_y_y[i] = 4.0 * g_xy_yy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_y_z[i] = 4.0 * g_xy_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_0_x_0_y_yy_z_x, g_x_0_x_0_y_yy_z_y, g_x_0_x_0_y_yy_z_z, g_xy_yy_xz_x, g_xy_yy_xz_y, g_xy_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_yy_z_x[i] = 4.0 * g_xy_yy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_z_y[i] = 4.0 * g_xy_yy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_z_z[i] = 4.0 * g_xy_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_0_x_0_y_yz_x_x, g_x_0_x_0_y_yz_x_y, g_x_0_x_0_y_yz_x_z, g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xy_yz_xx_x, g_xy_yz_xx_y, g_xy_yz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_yz_x_x[i] = -2.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xy_yz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_x_y[i] = -2.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xy_yz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_x_z[i] = -2.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xy_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_0_x_0_y_yz_y_x, g_x_0_x_0_y_yz_y_y, g_x_0_x_0_y_yz_y_z, g_xy_yz_xy_x, g_xy_yz_xy_y, g_xy_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_yz_y_x[i] = 4.0 * g_xy_yz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_y_y[i] = 4.0 * g_xy_yz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_y_z[i] = 4.0 * g_xy_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_0_x_0_y_yz_z_x, g_x_0_x_0_y_yz_z_y, g_x_0_x_0_y_yz_z_z, g_xy_yz_xz_x, g_xy_yz_xz_y, g_xy_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_yz_z_x[i] = 4.0 * g_xy_yz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_z_y[i] = 4.0 * g_xy_yz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_z_z[i] = 4.0 * g_xy_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_0_x_0_y_zz_x_x, g_x_0_x_0_y_zz_x_y, g_x_0_x_0_y_zz_x_z, g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xy_zz_xx_x, g_xy_zz_xx_y, g_xy_zz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_zz_x_x[i] = -2.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xy_zz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_x_y[i] = -2.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xy_zz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_x_z[i] = -2.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xy_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_0_x_0_y_zz_y_x, g_x_0_x_0_y_zz_y_y, g_x_0_x_0_y_zz_y_z, g_xy_zz_xy_x, g_xy_zz_xy_y, g_xy_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_zz_y_x[i] = 4.0 * g_xy_zz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_y_y[i] = 4.0 * g_xy_zz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_y_z[i] = 4.0 * g_xy_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_0_x_0_y_zz_z_x, g_x_0_x_0_y_zz_z_y, g_x_0_x_0_y_zz_z_z, g_xy_zz_xz_x, g_xy_zz_xz_y, g_xy_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_zz_z_x[i] = 4.0 * g_xy_zz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_z_y[i] = 4.0 * g_xy_zz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_z_z[i] = 4.0 * g_xy_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_0_x_0_z_xx_x_x, g_x_0_x_0_z_xx_x_y, g_x_0_x_0_z_xx_x_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_xz_xx_xx_x, g_xz_xx_xx_y, g_xz_xx_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xx_x_x[i] = -2.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xz_xx_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_x_y[i] = -2.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xz_xx_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_x_z[i] = -2.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xz_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_0_x_0_z_xx_y_x, g_x_0_x_0_z_xx_y_y, g_x_0_x_0_z_xx_y_z, g_xz_xx_xy_x, g_xz_xx_xy_y, g_xz_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xx_y_x[i] = 4.0 * g_xz_xx_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_y_y[i] = 4.0 * g_xz_xx_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_y_z[i] = 4.0 * g_xz_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_0_x_0_z_xx_z_x, g_x_0_x_0_z_xx_z_y, g_x_0_x_0_z_xx_z_z, g_xz_xx_xz_x, g_xz_xx_xz_y, g_xz_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xx_z_x[i] = 4.0 * g_xz_xx_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_z_y[i] = 4.0 * g_xz_xx_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_z_z[i] = 4.0 * g_xz_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_0_x_0_z_xy_x_x, g_x_0_x_0_z_xy_x_y, g_x_0_x_0_z_xy_x_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_xz_xy_xx_x, g_xz_xy_xx_y, g_xz_xy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xy_x_x[i] = -2.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xz_xy_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_x_y[i] = -2.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xz_xy_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_x_z[i] = -2.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xz_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_0_x_0_z_xy_y_x, g_x_0_x_0_z_xy_y_y, g_x_0_x_0_z_xy_y_z, g_xz_xy_xy_x, g_xz_xy_xy_y, g_xz_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xy_y_x[i] = 4.0 * g_xz_xy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_y_y[i] = 4.0 * g_xz_xy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_y_z[i] = 4.0 * g_xz_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_0_x_0_z_xy_z_x, g_x_0_x_0_z_xy_z_y, g_x_0_x_0_z_xy_z_z, g_xz_xy_xz_x, g_xz_xy_xz_y, g_xz_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xy_z_x[i] = 4.0 * g_xz_xy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_z_y[i] = 4.0 * g_xz_xy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_z_z[i] = 4.0 * g_xz_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_0_x_0_z_xz_x_x, g_x_0_x_0_z_xz_x_y, g_x_0_x_0_z_xz_x_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_xz_xz_xx_x, g_xz_xz_xx_y, g_xz_xz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xz_x_x[i] = -2.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xz_xz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_x_y[i] = -2.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xz_xz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_x_z[i] = -2.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xz_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_0_x_0_z_xz_y_x, g_x_0_x_0_z_xz_y_y, g_x_0_x_0_z_xz_y_z, g_xz_xz_xy_x, g_xz_xz_xy_y, g_xz_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xz_y_x[i] = 4.0 * g_xz_xz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_y_y[i] = 4.0 * g_xz_xz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_y_z[i] = 4.0 * g_xz_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_0_x_0_z_xz_z_x, g_x_0_x_0_z_xz_z_y, g_x_0_x_0_z_xz_z_z, g_xz_xz_xz_x, g_xz_xz_xz_y, g_xz_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xz_z_x[i] = 4.0 * g_xz_xz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_z_y[i] = 4.0 * g_xz_xz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_z_z[i] = 4.0 * g_xz_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_0_x_0_z_yy_x_x, g_x_0_x_0_z_yy_x_y, g_x_0_x_0_z_yy_x_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_xz_yy_xx_x, g_xz_yy_xx_y, g_xz_yy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_yy_x_x[i] = -2.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xz_yy_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_x_y[i] = -2.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xz_yy_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_x_z[i] = -2.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xz_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_0_x_0_z_yy_y_x, g_x_0_x_0_z_yy_y_y, g_x_0_x_0_z_yy_y_z, g_xz_yy_xy_x, g_xz_yy_xy_y, g_xz_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_yy_y_x[i] = 4.0 * g_xz_yy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_y_y[i] = 4.0 * g_xz_yy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_y_z[i] = 4.0 * g_xz_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_0_x_0_z_yy_z_x, g_x_0_x_0_z_yy_z_y, g_x_0_x_0_z_yy_z_z, g_xz_yy_xz_x, g_xz_yy_xz_y, g_xz_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_yy_z_x[i] = 4.0 * g_xz_yy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_z_y[i] = 4.0 * g_xz_yy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_z_z[i] = 4.0 * g_xz_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_0_x_0_z_yz_x_x, g_x_0_x_0_z_yz_x_y, g_x_0_x_0_z_yz_x_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_xz_yz_xx_x, g_xz_yz_xx_y, g_xz_yz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_yz_x_x[i] = -2.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xz_yz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_x_y[i] = -2.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xz_yz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_x_z[i] = -2.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xz_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_0_x_0_z_yz_y_x, g_x_0_x_0_z_yz_y_y, g_x_0_x_0_z_yz_y_z, g_xz_yz_xy_x, g_xz_yz_xy_y, g_xz_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_yz_y_x[i] = 4.0 * g_xz_yz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_y_y[i] = 4.0 * g_xz_yz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_y_z[i] = 4.0 * g_xz_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_0_x_0_z_yz_z_x, g_x_0_x_0_z_yz_z_y, g_x_0_x_0_z_yz_z_z, g_xz_yz_xz_x, g_xz_yz_xz_y, g_xz_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_yz_z_x[i] = 4.0 * g_xz_yz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_z_y[i] = 4.0 * g_xz_yz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_z_z[i] = 4.0 * g_xz_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_0_x_0_z_zz_x_x, g_x_0_x_0_z_zz_x_y, g_x_0_x_0_z_zz_x_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_xz_zz_xx_x, g_xz_zz_xx_y, g_xz_zz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_zz_x_x[i] = -2.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xz_zz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_x_y[i] = -2.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xz_zz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_x_z[i] = -2.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xz_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_0_x_0_z_zz_y_x, g_x_0_x_0_z_zz_y_y, g_x_0_x_0_z_zz_y_z, g_xz_zz_xy_x, g_xz_zz_xy_y, g_xz_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_zz_y_x[i] = 4.0 * g_xz_zz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_y_y[i] = 4.0 * g_xz_zz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_y_z[i] = 4.0 * g_xz_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_0_x_0_z_zz_z_x, g_x_0_x_0_z_zz_z_y, g_x_0_x_0_z_zz_z_z, g_xz_zz_xz_x, g_xz_zz_xz_y, g_xz_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_zz_z_x[i] = 4.0 * g_xz_zz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_z_y[i] = 4.0 * g_xz_zz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_z_z[i] = 4.0 * g_xz_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_0_xx_xy_x, g_0_xx_xy_y, g_0_xx_xy_z, g_x_0_y_0_x_xx_x_x, g_x_0_y_0_x_xx_x_y, g_x_0_y_0_x_xx_x_z, g_xx_xx_xy_x, g_xx_xx_xy_y, g_xx_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xx_x_x[i] = -2.0 * g_0_xx_xy_x[i] * c_exps[i] + 4.0 * g_xx_xx_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_x_y[i] = -2.0 * g_0_xx_xy_y[i] * c_exps[i] + 4.0 * g_xx_xx_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_x_z[i] = -2.0 * g_0_xx_xy_z[i] * c_exps[i] + 4.0 * g_xx_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_0_xx_yy_x, g_0_xx_yy_y, g_0_xx_yy_z, g_x_0_y_0_x_xx_y_x, g_x_0_y_0_x_xx_y_y, g_x_0_y_0_x_xx_y_z, g_xx_xx_0_x, g_xx_xx_0_y, g_xx_xx_0_z, g_xx_xx_yy_x, g_xx_xx_yy_y, g_xx_xx_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xx_y_x[i] = g_0_xx_0_x[i] - 2.0 * g_0_xx_yy_x[i] * c_exps[i] - 2.0 * g_xx_xx_0_x[i] * a_exp + 4.0 * g_xx_xx_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_y_y[i] = g_0_xx_0_y[i] - 2.0 * g_0_xx_yy_y[i] * c_exps[i] - 2.0 * g_xx_xx_0_y[i] * a_exp + 4.0 * g_xx_xx_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_y_z[i] = g_0_xx_0_z[i] - 2.0 * g_0_xx_yy_z[i] * c_exps[i] - 2.0 * g_xx_xx_0_z[i] * a_exp + 4.0 * g_xx_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_0_xx_yz_x, g_0_xx_yz_y, g_0_xx_yz_z, g_x_0_y_0_x_xx_z_x, g_x_0_y_0_x_xx_z_y, g_x_0_y_0_x_xx_z_z, g_xx_xx_yz_x, g_xx_xx_yz_y, g_xx_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xx_z_x[i] = -2.0 * g_0_xx_yz_x[i] * c_exps[i] + 4.0 * g_xx_xx_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_z_y[i] = -2.0 * g_0_xx_yz_y[i] * c_exps[i] + 4.0 * g_xx_xx_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_z_z[i] = -2.0 * g_0_xx_yz_z[i] * c_exps[i] + 4.0 * g_xx_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_0_xy_xy_x, g_0_xy_xy_y, g_0_xy_xy_z, g_x_0_y_0_x_xy_x_x, g_x_0_y_0_x_xy_x_y, g_x_0_y_0_x_xy_x_z, g_xx_xy_xy_x, g_xx_xy_xy_y, g_xx_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xy_x_x[i] = -2.0 * g_0_xy_xy_x[i] * c_exps[i] + 4.0 * g_xx_xy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_x_y[i] = -2.0 * g_0_xy_xy_y[i] * c_exps[i] + 4.0 * g_xx_xy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_x_z[i] = -2.0 * g_0_xy_xy_z[i] * c_exps[i] + 4.0 * g_xx_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_0_xy_yy_x, g_0_xy_yy_y, g_0_xy_yy_z, g_x_0_y_0_x_xy_y_x, g_x_0_y_0_x_xy_y_y, g_x_0_y_0_x_xy_y_z, g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z, g_xx_xy_yy_x, g_xx_xy_yy_y, g_xx_xy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xy_y_x[i] = g_0_xy_0_x[i] - 2.0 * g_0_xy_yy_x[i] * c_exps[i] - 2.0 * g_xx_xy_0_x[i] * a_exp + 4.0 * g_xx_xy_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_y_y[i] = g_0_xy_0_y[i] - 2.0 * g_0_xy_yy_y[i] * c_exps[i] - 2.0 * g_xx_xy_0_y[i] * a_exp + 4.0 * g_xx_xy_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_y_z[i] = g_0_xy_0_z[i] - 2.0 * g_0_xy_yy_z[i] * c_exps[i] - 2.0 * g_xx_xy_0_z[i] * a_exp + 4.0 * g_xx_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_0_xy_yz_x, g_0_xy_yz_y, g_0_xy_yz_z, g_x_0_y_0_x_xy_z_x, g_x_0_y_0_x_xy_z_y, g_x_0_y_0_x_xy_z_z, g_xx_xy_yz_x, g_xx_xy_yz_y, g_xx_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xy_z_x[i] = -2.0 * g_0_xy_yz_x[i] * c_exps[i] + 4.0 * g_xx_xy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_z_y[i] = -2.0 * g_0_xy_yz_y[i] * c_exps[i] + 4.0 * g_xx_xy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_z_z[i] = -2.0 * g_0_xy_yz_z[i] * c_exps[i] + 4.0 * g_xx_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_0_xz_xy_x, g_0_xz_xy_y, g_0_xz_xy_z, g_x_0_y_0_x_xz_x_x, g_x_0_y_0_x_xz_x_y, g_x_0_y_0_x_xz_x_z, g_xx_xz_xy_x, g_xx_xz_xy_y, g_xx_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xz_x_x[i] = -2.0 * g_0_xz_xy_x[i] * c_exps[i] + 4.0 * g_xx_xz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_x_y[i] = -2.0 * g_0_xz_xy_y[i] * c_exps[i] + 4.0 * g_xx_xz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_x_z[i] = -2.0 * g_0_xz_xy_z[i] * c_exps[i] + 4.0 * g_xx_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_0_xz_yy_x, g_0_xz_yy_y, g_0_xz_yy_z, g_x_0_y_0_x_xz_y_x, g_x_0_y_0_x_xz_y_y, g_x_0_y_0_x_xz_y_z, g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z, g_xx_xz_yy_x, g_xx_xz_yy_y, g_xx_xz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xz_y_x[i] = g_0_xz_0_x[i] - 2.0 * g_0_xz_yy_x[i] * c_exps[i] - 2.0 * g_xx_xz_0_x[i] * a_exp + 4.0 * g_xx_xz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_y_y[i] = g_0_xz_0_y[i] - 2.0 * g_0_xz_yy_y[i] * c_exps[i] - 2.0 * g_xx_xz_0_y[i] * a_exp + 4.0 * g_xx_xz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_y_z[i] = g_0_xz_0_z[i] - 2.0 * g_0_xz_yy_z[i] * c_exps[i] - 2.0 * g_xx_xz_0_z[i] * a_exp + 4.0 * g_xx_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_0_xz_yz_x, g_0_xz_yz_y, g_0_xz_yz_z, g_x_0_y_0_x_xz_z_x, g_x_0_y_0_x_xz_z_y, g_x_0_y_0_x_xz_z_z, g_xx_xz_yz_x, g_xx_xz_yz_y, g_xx_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xz_z_x[i] = -2.0 * g_0_xz_yz_x[i] * c_exps[i] + 4.0 * g_xx_xz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_z_y[i] = -2.0 * g_0_xz_yz_y[i] * c_exps[i] + 4.0 * g_xx_xz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_z_z[i] = -2.0 * g_0_xz_yz_z[i] * c_exps[i] + 4.0 * g_xx_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_0_yy_xy_x, g_0_yy_xy_y, g_0_yy_xy_z, g_x_0_y_0_x_yy_x_x, g_x_0_y_0_x_yy_x_y, g_x_0_y_0_x_yy_x_z, g_xx_yy_xy_x, g_xx_yy_xy_y, g_xx_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_yy_x_x[i] = -2.0 * g_0_yy_xy_x[i] * c_exps[i] + 4.0 * g_xx_yy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_x_y[i] = -2.0 * g_0_yy_xy_y[i] * c_exps[i] + 4.0 * g_xx_yy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_x_z[i] = -2.0 * g_0_yy_xy_z[i] * c_exps[i] + 4.0 * g_xx_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_0_yy_yy_x, g_0_yy_yy_y, g_0_yy_yy_z, g_x_0_y_0_x_yy_y_x, g_x_0_y_0_x_yy_y_y, g_x_0_y_0_x_yy_y_z, g_xx_yy_0_x, g_xx_yy_0_y, g_xx_yy_0_z, g_xx_yy_yy_x, g_xx_yy_yy_y, g_xx_yy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_yy_y_x[i] = g_0_yy_0_x[i] - 2.0 * g_0_yy_yy_x[i] * c_exps[i] - 2.0 * g_xx_yy_0_x[i] * a_exp + 4.0 * g_xx_yy_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_y_y[i] = g_0_yy_0_y[i] - 2.0 * g_0_yy_yy_y[i] * c_exps[i] - 2.0 * g_xx_yy_0_y[i] * a_exp + 4.0 * g_xx_yy_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_y_z[i] = g_0_yy_0_z[i] - 2.0 * g_0_yy_yy_z[i] * c_exps[i] - 2.0 * g_xx_yy_0_z[i] * a_exp + 4.0 * g_xx_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_0_yy_yz_x, g_0_yy_yz_y, g_0_yy_yz_z, g_x_0_y_0_x_yy_z_x, g_x_0_y_0_x_yy_z_y, g_x_0_y_0_x_yy_z_z, g_xx_yy_yz_x, g_xx_yy_yz_y, g_xx_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_yy_z_x[i] = -2.0 * g_0_yy_yz_x[i] * c_exps[i] + 4.0 * g_xx_yy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_z_y[i] = -2.0 * g_0_yy_yz_y[i] * c_exps[i] + 4.0 * g_xx_yy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_z_z[i] = -2.0 * g_0_yy_yz_z[i] * c_exps[i] + 4.0 * g_xx_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_0_yz_xy_x, g_0_yz_xy_y, g_0_yz_xy_z, g_x_0_y_0_x_yz_x_x, g_x_0_y_0_x_yz_x_y, g_x_0_y_0_x_yz_x_z, g_xx_yz_xy_x, g_xx_yz_xy_y, g_xx_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_yz_x_x[i] = -2.0 * g_0_yz_xy_x[i] * c_exps[i] + 4.0 * g_xx_yz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_x_y[i] = -2.0 * g_0_yz_xy_y[i] * c_exps[i] + 4.0 * g_xx_yz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_x_z[i] = -2.0 * g_0_yz_xy_z[i] * c_exps[i] + 4.0 * g_xx_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_0_yz_yy_x, g_0_yz_yy_y, g_0_yz_yy_z, g_x_0_y_0_x_yz_y_x, g_x_0_y_0_x_yz_y_y, g_x_0_y_0_x_yz_y_z, g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z, g_xx_yz_yy_x, g_xx_yz_yy_y, g_xx_yz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_yz_y_x[i] = g_0_yz_0_x[i] - 2.0 * g_0_yz_yy_x[i] * c_exps[i] - 2.0 * g_xx_yz_0_x[i] * a_exp + 4.0 * g_xx_yz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_y_y[i] = g_0_yz_0_y[i] - 2.0 * g_0_yz_yy_y[i] * c_exps[i] - 2.0 * g_xx_yz_0_y[i] * a_exp + 4.0 * g_xx_yz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_y_z[i] = g_0_yz_0_z[i] - 2.0 * g_0_yz_yy_z[i] * c_exps[i] - 2.0 * g_xx_yz_0_z[i] * a_exp + 4.0 * g_xx_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_0_yz_yz_x, g_0_yz_yz_y, g_0_yz_yz_z, g_x_0_y_0_x_yz_z_x, g_x_0_y_0_x_yz_z_y, g_x_0_y_0_x_yz_z_z, g_xx_yz_yz_x, g_xx_yz_yz_y, g_xx_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_yz_z_x[i] = -2.0 * g_0_yz_yz_x[i] * c_exps[i] + 4.0 * g_xx_yz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_z_y[i] = -2.0 * g_0_yz_yz_y[i] * c_exps[i] + 4.0 * g_xx_yz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_z_z[i] = -2.0 * g_0_yz_yz_z[i] * c_exps[i] + 4.0 * g_xx_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_0_zz_xy_x, g_0_zz_xy_y, g_0_zz_xy_z, g_x_0_y_0_x_zz_x_x, g_x_0_y_0_x_zz_x_y, g_x_0_y_0_x_zz_x_z, g_xx_zz_xy_x, g_xx_zz_xy_y, g_xx_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_zz_x_x[i] = -2.0 * g_0_zz_xy_x[i] * c_exps[i] + 4.0 * g_xx_zz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_x_y[i] = -2.0 * g_0_zz_xy_y[i] * c_exps[i] + 4.0 * g_xx_zz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_x_z[i] = -2.0 * g_0_zz_xy_z[i] * c_exps[i] + 4.0 * g_xx_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_0_zz_yy_x, g_0_zz_yy_y, g_0_zz_yy_z, g_x_0_y_0_x_zz_y_x, g_x_0_y_0_x_zz_y_y, g_x_0_y_0_x_zz_y_z, g_xx_zz_0_x, g_xx_zz_0_y, g_xx_zz_0_z, g_xx_zz_yy_x, g_xx_zz_yy_y, g_xx_zz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_zz_y_x[i] = g_0_zz_0_x[i] - 2.0 * g_0_zz_yy_x[i] * c_exps[i] - 2.0 * g_xx_zz_0_x[i] * a_exp + 4.0 * g_xx_zz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_y_y[i] = g_0_zz_0_y[i] - 2.0 * g_0_zz_yy_y[i] * c_exps[i] - 2.0 * g_xx_zz_0_y[i] * a_exp + 4.0 * g_xx_zz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_y_z[i] = g_0_zz_0_z[i] - 2.0 * g_0_zz_yy_z[i] * c_exps[i] - 2.0 * g_xx_zz_0_z[i] * a_exp + 4.0 * g_xx_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_0_zz_yz_x, g_0_zz_yz_y, g_0_zz_yz_z, g_x_0_y_0_x_zz_z_x, g_x_0_y_0_x_zz_z_y, g_x_0_y_0_x_zz_z_z, g_xx_zz_yz_x, g_xx_zz_yz_y, g_xx_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_zz_z_x[i] = -2.0 * g_0_zz_yz_x[i] * c_exps[i] + 4.0 * g_xx_zz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_z_y[i] = -2.0 * g_0_zz_yz_y[i] * c_exps[i] + 4.0 * g_xx_zz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_z_z[i] = -2.0 * g_0_zz_yz_z[i] * c_exps[i] + 4.0 * g_xx_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_x_0_y_0_y_xx_x_x, g_x_0_y_0_y_xx_x_y, g_x_0_y_0_y_xx_x_z, g_xy_xx_xy_x, g_xy_xx_xy_y, g_xy_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xx_x_x[i] = 4.0 * g_xy_xx_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_x_y[i] = 4.0 * g_xy_xx_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_x_z[i] = 4.0 * g_xy_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_x_0_y_0_y_xx_y_x, g_x_0_y_0_y_xx_y_y, g_x_0_y_0_y_xx_y_z, g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xy_xx_yy_x, g_xy_xx_yy_y, g_xy_xx_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xx_y_x[i] = -2.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xy_xx_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_y_y[i] = -2.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xy_xx_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_y_z[i] = -2.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xy_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_x_0_y_0_y_xx_z_x, g_x_0_y_0_y_xx_z_y, g_x_0_y_0_y_xx_z_z, g_xy_xx_yz_x, g_xy_xx_yz_y, g_xy_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xx_z_x[i] = 4.0 * g_xy_xx_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_z_y[i] = 4.0 * g_xy_xx_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_z_z[i] = 4.0 * g_xy_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_x_0_y_0_y_xy_x_x, g_x_0_y_0_y_xy_x_y, g_x_0_y_0_y_xy_x_z, g_xy_xy_xy_x, g_xy_xy_xy_y, g_xy_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xy_x_x[i] = 4.0 * g_xy_xy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_x_y[i] = 4.0 * g_xy_xy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_x_z[i] = 4.0 * g_xy_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_x_0_y_0_y_xy_y_x, g_x_0_y_0_y_xy_y_y, g_x_0_y_0_y_xy_y_z, g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xy_xy_yy_x, g_xy_xy_yy_y, g_xy_xy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xy_y_x[i] = -2.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xy_xy_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_y_y[i] = -2.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xy_xy_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_y_z[i] = -2.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xy_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_x_0_y_0_y_xy_z_x, g_x_0_y_0_y_xy_z_y, g_x_0_y_0_y_xy_z_z, g_xy_xy_yz_x, g_xy_xy_yz_y, g_xy_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xy_z_x[i] = 4.0 * g_xy_xy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_z_y[i] = 4.0 * g_xy_xy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_z_z[i] = 4.0 * g_xy_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_x_0_y_0_y_xz_x_x, g_x_0_y_0_y_xz_x_y, g_x_0_y_0_y_xz_x_z, g_xy_xz_xy_x, g_xy_xz_xy_y, g_xy_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xz_x_x[i] = 4.0 * g_xy_xz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_x_y[i] = 4.0 * g_xy_xz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_x_z[i] = 4.0 * g_xy_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_x_0_y_0_y_xz_y_x, g_x_0_y_0_y_xz_y_y, g_x_0_y_0_y_xz_y_z, g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xy_xz_yy_x, g_xy_xz_yy_y, g_xy_xz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xz_y_x[i] = -2.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xy_xz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_y_y[i] = -2.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xy_xz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_y_z[i] = -2.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xy_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_x_0_y_0_y_xz_z_x, g_x_0_y_0_y_xz_z_y, g_x_0_y_0_y_xz_z_z, g_xy_xz_yz_x, g_xy_xz_yz_y, g_xy_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xz_z_x[i] = 4.0 * g_xy_xz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_z_y[i] = 4.0 * g_xy_xz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_z_z[i] = 4.0 * g_xy_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_x_0_y_0_y_yy_x_x, g_x_0_y_0_y_yy_x_y, g_x_0_y_0_y_yy_x_z, g_xy_yy_xy_x, g_xy_yy_xy_y, g_xy_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_yy_x_x[i] = 4.0 * g_xy_yy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_x_y[i] = 4.0 * g_xy_yy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_x_z[i] = 4.0 * g_xy_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_x_0_y_0_y_yy_y_x, g_x_0_y_0_y_yy_y_y, g_x_0_y_0_y_yy_y_z, g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xy_yy_yy_x, g_xy_yy_yy_y, g_xy_yy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_yy_y_x[i] = -2.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xy_yy_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_y_y[i] = -2.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xy_yy_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_y_z[i] = -2.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xy_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_x_0_y_0_y_yy_z_x, g_x_0_y_0_y_yy_z_y, g_x_0_y_0_y_yy_z_z, g_xy_yy_yz_x, g_xy_yy_yz_y, g_xy_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_yy_z_x[i] = 4.0 * g_xy_yy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_z_y[i] = 4.0 * g_xy_yy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_z_z[i] = 4.0 * g_xy_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_x_0_y_0_y_yz_x_x, g_x_0_y_0_y_yz_x_y, g_x_0_y_0_y_yz_x_z, g_xy_yz_xy_x, g_xy_yz_xy_y, g_xy_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_yz_x_x[i] = 4.0 * g_xy_yz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_x_y[i] = 4.0 * g_xy_yz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_x_z[i] = 4.0 * g_xy_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_x_0_y_0_y_yz_y_x, g_x_0_y_0_y_yz_y_y, g_x_0_y_0_y_yz_y_z, g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xy_yz_yy_x, g_xy_yz_yy_y, g_xy_yz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_yz_y_x[i] = -2.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xy_yz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_y_y[i] = -2.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xy_yz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_y_z[i] = -2.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xy_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_x_0_y_0_y_yz_z_x, g_x_0_y_0_y_yz_z_y, g_x_0_y_0_y_yz_z_z, g_xy_yz_yz_x, g_xy_yz_yz_y, g_xy_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_yz_z_x[i] = 4.0 * g_xy_yz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_z_y[i] = 4.0 * g_xy_yz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_z_z[i] = 4.0 * g_xy_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_x_0_y_0_y_zz_x_x, g_x_0_y_0_y_zz_x_y, g_x_0_y_0_y_zz_x_z, g_xy_zz_xy_x, g_xy_zz_xy_y, g_xy_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_zz_x_x[i] = 4.0 * g_xy_zz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_x_y[i] = 4.0 * g_xy_zz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_x_z[i] = 4.0 * g_xy_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_x_0_y_0_y_zz_y_x, g_x_0_y_0_y_zz_y_y, g_x_0_y_0_y_zz_y_z, g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xy_zz_yy_x, g_xy_zz_yy_y, g_xy_zz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_zz_y_x[i] = -2.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xy_zz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_y_y[i] = -2.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xy_zz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_y_z[i] = -2.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xy_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_x_0_y_0_y_zz_z_x, g_x_0_y_0_y_zz_z_y, g_x_0_y_0_y_zz_z_z, g_xy_zz_yz_x, g_xy_zz_yz_y, g_xy_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_zz_z_x[i] = 4.0 * g_xy_zz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_z_y[i] = 4.0 * g_xy_zz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_z_z[i] = 4.0 * g_xy_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_x_0_y_0_z_xx_x_x, g_x_0_y_0_z_xx_x_y, g_x_0_y_0_z_xx_x_z, g_xz_xx_xy_x, g_xz_xx_xy_y, g_xz_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xx_x_x[i] = 4.0 * g_xz_xx_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_x_y[i] = 4.0 * g_xz_xx_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_x_z[i] = 4.0 * g_xz_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_x_0_y_0_z_xx_y_x, g_x_0_y_0_z_xx_y_y, g_x_0_y_0_z_xx_y_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_xz_xx_yy_x, g_xz_xx_yy_y, g_xz_xx_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xx_y_x[i] = -2.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xz_xx_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_y_y[i] = -2.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xz_xx_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_y_z[i] = -2.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xz_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_x_0_y_0_z_xx_z_x, g_x_0_y_0_z_xx_z_y, g_x_0_y_0_z_xx_z_z, g_xz_xx_yz_x, g_xz_xx_yz_y, g_xz_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xx_z_x[i] = 4.0 * g_xz_xx_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_z_y[i] = 4.0 * g_xz_xx_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_z_z[i] = 4.0 * g_xz_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_x_0_y_0_z_xy_x_x, g_x_0_y_0_z_xy_x_y, g_x_0_y_0_z_xy_x_z, g_xz_xy_xy_x, g_xz_xy_xy_y, g_xz_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xy_x_x[i] = 4.0 * g_xz_xy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_x_y[i] = 4.0 * g_xz_xy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_x_z[i] = 4.0 * g_xz_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_x_0_y_0_z_xy_y_x, g_x_0_y_0_z_xy_y_y, g_x_0_y_0_z_xy_y_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_xz_xy_yy_x, g_xz_xy_yy_y, g_xz_xy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xy_y_x[i] = -2.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xz_xy_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_y_y[i] = -2.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xz_xy_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_y_z[i] = -2.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xz_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_x_0_y_0_z_xy_z_x, g_x_0_y_0_z_xy_z_y, g_x_0_y_0_z_xy_z_z, g_xz_xy_yz_x, g_xz_xy_yz_y, g_xz_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xy_z_x[i] = 4.0 * g_xz_xy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_z_y[i] = 4.0 * g_xz_xy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_z_z[i] = 4.0 * g_xz_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_x_0_y_0_z_xz_x_x, g_x_0_y_0_z_xz_x_y, g_x_0_y_0_z_xz_x_z, g_xz_xz_xy_x, g_xz_xz_xy_y, g_xz_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xz_x_x[i] = 4.0 * g_xz_xz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_x_y[i] = 4.0 * g_xz_xz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_x_z[i] = 4.0 * g_xz_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_x_0_y_0_z_xz_y_x, g_x_0_y_0_z_xz_y_y, g_x_0_y_0_z_xz_y_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_xz_xz_yy_x, g_xz_xz_yy_y, g_xz_xz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xz_y_x[i] = -2.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xz_xz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_y_y[i] = -2.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xz_xz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_y_z[i] = -2.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xz_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_x_0_y_0_z_xz_z_x, g_x_0_y_0_z_xz_z_y, g_x_0_y_0_z_xz_z_z, g_xz_xz_yz_x, g_xz_xz_yz_y, g_xz_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xz_z_x[i] = 4.0 * g_xz_xz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_z_y[i] = 4.0 * g_xz_xz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_z_z[i] = 4.0 * g_xz_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_x_0_y_0_z_yy_x_x, g_x_0_y_0_z_yy_x_y, g_x_0_y_0_z_yy_x_z, g_xz_yy_xy_x, g_xz_yy_xy_y, g_xz_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_yy_x_x[i] = 4.0 * g_xz_yy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_x_y[i] = 4.0 * g_xz_yy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_x_z[i] = 4.0 * g_xz_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_x_0_y_0_z_yy_y_x, g_x_0_y_0_z_yy_y_y, g_x_0_y_0_z_yy_y_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_xz_yy_yy_x, g_xz_yy_yy_y, g_xz_yy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_yy_y_x[i] = -2.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xz_yy_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_y_y[i] = -2.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xz_yy_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_y_z[i] = -2.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xz_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_x_0_y_0_z_yy_z_x, g_x_0_y_0_z_yy_z_y, g_x_0_y_0_z_yy_z_z, g_xz_yy_yz_x, g_xz_yy_yz_y, g_xz_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_yy_z_x[i] = 4.0 * g_xz_yy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_z_y[i] = 4.0 * g_xz_yy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_z_z[i] = 4.0 * g_xz_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_x_0_y_0_z_yz_x_x, g_x_0_y_0_z_yz_x_y, g_x_0_y_0_z_yz_x_z, g_xz_yz_xy_x, g_xz_yz_xy_y, g_xz_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_yz_x_x[i] = 4.0 * g_xz_yz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_x_y[i] = 4.0 * g_xz_yz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_x_z[i] = 4.0 * g_xz_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_x_0_y_0_z_yz_y_x, g_x_0_y_0_z_yz_y_y, g_x_0_y_0_z_yz_y_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_xz_yz_yy_x, g_xz_yz_yy_y, g_xz_yz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_yz_y_x[i] = -2.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xz_yz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_y_y[i] = -2.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xz_yz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_y_z[i] = -2.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xz_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_x_0_y_0_z_yz_z_x, g_x_0_y_0_z_yz_z_y, g_x_0_y_0_z_yz_z_z, g_xz_yz_yz_x, g_xz_yz_yz_y, g_xz_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_yz_z_x[i] = 4.0 * g_xz_yz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_z_y[i] = 4.0 * g_xz_yz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_z_z[i] = 4.0 * g_xz_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_x_0_y_0_z_zz_x_x, g_x_0_y_0_z_zz_x_y, g_x_0_y_0_z_zz_x_z, g_xz_zz_xy_x, g_xz_zz_xy_y, g_xz_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_zz_x_x[i] = 4.0 * g_xz_zz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_x_y[i] = 4.0 * g_xz_zz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_x_z[i] = 4.0 * g_xz_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_x_0_y_0_z_zz_y_x, g_x_0_y_0_z_zz_y_y, g_x_0_y_0_z_zz_y_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_xz_zz_yy_x, g_xz_zz_yy_y, g_xz_zz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_zz_y_x[i] = -2.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xz_zz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_y_y[i] = -2.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xz_zz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_y_z[i] = -2.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xz_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_x_0_y_0_z_zz_z_x, g_x_0_y_0_z_zz_z_y, g_x_0_y_0_z_zz_z_z, g_xz_zz_yz_x, g_xz_zz_yz_y, g_xz_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_zz_z_x[i] = 4.0 * g_xz_zz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_z_y[i] = 4.0 * g_xz_zz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_z_z[i] = 4.0 * g_xz_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_0_xx_xz_x, g_0_xx_xz_y, g_0_xx_xz_z, g_x_0_z_0_x_xx_x_x, g_x_0_z_0_x_xx_x_y, g_x_0_z_0_x_xx_x_z, g_xx_xx_xz_x, g_xx_xx_xz_y, g_xx_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xx_x_x[i] = -2.0 * g_0_xx_xz_x[i] * c_exps[i] + 4.0 * g_xx_xx_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_x_y[i] = -2.0 * g_0_xx_xz_y[i] * c_exps[i] + 4.0 * g_xx_xx_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_x_z[i] = -2.0 * g_0_xx_xz_z[i] * c_exps[i] + 4.0 * g_xx_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_0_xx_yz_x, g_0_xx_yz_y, g_0_xx_yz_z, g_x_0_z_0_x_xx_y_x, g_x_0_z_0_x_xx_y_y, g_x_0_z_0_x_xx_y_z, g_xx_xx_yz_x, g_xx_xx_yz_y, g_xx_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xx_y_x[i] = -2.0 * g_0_xx_yz_x[i] * c_exps[i] + 4.0 * g_xx_xx_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_y_y[i] = -2.0 * g_0_xx_yz_y[i] * c_exps[i] + 4.0 * g_xx_xx_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_y_z[i] = -2.0 * g_0_xx_yz_z[i] * c_exps[i] + 4.0 * g_xx_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_0_xx_zz_x, g_0_xx_zz_y, g_0_xx_zz_z, g_x_0_z_0_x_xx_z_x, g_x_0_z_0_x_xx_z_y, g_x_0_z_0_x_xx_z_z, g_xx_xx_0_x, g_xx_xx_0_y, g_xx_xx_0_z, g_xx_xx_zz_x, g_xx_xx_zz_y, g_xx_xx_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xx_z_x[i] = g_0_xx_0_x[i] - 2.0 * g_0_xx_zz_x[i] * c_exps[i] - 2.0 * g_xx_xx_0_x[i] * a_exp + 4.0 * g_xx_xx_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_z_y[i] = g_0_xx_0_y[i] - 2.0 * g_0_xx_zz_y[i] * c_exps[i] - 2.0 * g_xx_xx_0_y[i] * a_exp + 4.0 * g_xx_xx_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_z_z[i] = g_0_xx_0_z[i] - 2.0 * g_0_xx_zz_z[i] * c_exps[i] - 2.0 * g_xx_xx_0_z[i] * a_exp + 4.0 * g_xx_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_0_xy_xz_x, g_0_xy_xz_y, g_0_xy_xz_z, g_x_0_z_0_x_xy_x_x, g_x_0_z_0_x_xy_x_y, g_x_0_z_0_x_xy_x_z, g_xx_xy_xz_x, g_xx_xy_xz_y, g_xx_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xy_x_x[i] = -2.0 * g_0_xy_xz_x[i] * c_exps[i] + 4.0 * g_xx_xy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_x_y[i] = -2.0 * g_0_xy_xz_y[i] * c_exps[i] + 4.0 * g_xx_xy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_x_z[i] = -2.0 * g_0_xy_xz_z[i] * c_exps[i] + 4.0 * g_xx_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_0_xy_yz_x, g_0_xy_yz_y, g_0_xy_yz_z, g_x_0_z_0_x_xy_y_x, g_x_0_z_0_x_xy_y_y, g_x_0_z_0_x_xy_y_z, g_xx_xy_yz_x, g_xx_xy_yz_y, g_xx_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xy_y_x[i] = -2.0 * g_0_xy_yz_x[i] * c_exps[i] + 4.0 * g_xx_xy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_y_y[i] = -2.0 * g_0_xy_yz_y[i] * c_exps[i] + 4.0 * g_xx_xy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_y_z[i] = -2.0 * g_0_xy_yz_z[i] * c_exps[i] + 4.0 * g_xx_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_0_xy_zz_x, g_0_xy_zz_y, g_0_xy_zz_z, g_x_0_z_0_x_xy_z_x, g_x_0_z_0_x_xy_z_y, g_x_0_z_0_x_xy_z_z, g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z, g_xx_xy_zz_x, g_xx_xy_zz_y, g_xx_xy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xy_z_x[i] = g_0_xy_0_x[i] - 2.0 * g_0_xy_zz_x[i] * c_exps[i] - 2.0 * g_xx_xy_0_x[i] * a_exp + 4.0 * g_xx_xy_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_z_y[i] = g_0_xy_0_y[i] - 2.0 * g_0_xy_zz_y[i] * c_exps[i] - 2.0 * g_xx_xy_0_y[i] * a_exp + 4.0 * g_xx_xy_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_z_z[i] = g_0_xy_0_z[i] - 2.0 * g_0_xy_zz_z[i] * c_exps[i] - 2.0 * g_xx_xy_0_z[i] * a_exp + 4.0 * g_xx_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_0_xz_xz_x, g_0_xz_xz_y, g_0_xz_xz_z, g_x_0_z_0_x_xz_x_x, g_x_0_z_0_x_xz_x_y, g_x_0_z_0_x_xz_x_z, g_xx_xz_xz_x, g_xx_xz_xz_y, g_xx_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xz_x_x[i] = -2.0 * g_0_xz_xz_x[i] * c_exps[i] + 4.0 * g_xx_xz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_x_y[i] = -2.0 * g_0_xz_xz_y[i] * c_exps[i] + 4.0 * g_xx_xz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_x_z[i] = -2.0 * g_0_xz_xz_z[i] * c_exps[i] + 4.0 * g_xx_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_0_xz_yz_x, g_0_xz_yz_y, g_0_xz_yz_z, g_x_0_z_0_x_xz_y_x, g_x_0_z_0_x_xz_y_y, g_x_0_z_0_x_xz_y_z, g_xx_xz_yz_x, g_xx_xz_yz_y, g_xx_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xz_y_x[i] = -2.0 * g_0_xz_yz_x[i] * c_exps[i] + 4.0 * g_xx_xz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_y_y[i] = -2.0 * g_0_xz_yz_y[i] * c_exps[i] + 4.0 * g_xx_xz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_y_z[i] = -2.0 * g_0_xz_yz_z[i] * c_exps[i] + 4.0 * g_xx_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_0_xz_zz_x, g_0_xz_zz_y, g_0_xz_zz_z, g_x_0_z_0_x_xz_z_x, g_x_0_z_0_x_xz_z_y, g_x_0_z_0_x_xz_z_z, g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z, g_xx_xz_zz_x, g_xx_xz_zz_y, g_xx_xz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xz_z_x[i] = g_0_xz_0_x[i] - 2.0 * g_0_xz_zz_x[i] * c_exps[i] - 2.0 * g_xx_xz_0_x[i] * a_exp + 4.0 * g_xx_xz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_z_y[i] = g_0_xz_0_y[i] - 2.0 * g_0_xz_zz_y[i] * c_exps[i] - 2.0 * g_xx_xz_0_y[i] * a_exp + 4.0 * g_xx_xz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_z_z[i] = g_0_xz_0_z[i] - 2.0 * g_0_xz_zz_z[i] * c_exps[i] - 2.0 * g_xx_xz_0_z[i] * a_exp + 4.0 * g_xx_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_0_yy_xz_x, g_0_yy_xz_y, g_0_yy_xz_z, g_x_0_z_0_x_yy_x_x, g_x_0_z_0_x_yy_x_y, g_x_0_z_0_x_yy_x_z, g_xx_yy_xz_x, g_xx_yy_xz_y, g_xx_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_yy_x_x[i] = -2.0 * g_0_yy_xz_x[i] * c_exps[i] + 4.0 * g_xx_yy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_x_y[i] = -2.0 * g_0_yy_xz_y[i] * c_exps[i] + 4.0 * g_xx_yy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_x_z[i] = -2.0 * g_0_yy_xz_z[i] * c_exps[i] + 4.0 * g_xx_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_0_yy_yz_x, g_0_yy_yz_y, g_0_yy_yz_z, g_x_0_z_0_x_yy_y_x, g_x_0_z_0_x_yy_y_y, g_x_0_z_0_x_yy_y_z, g_xx_yy_yz_x, g_xx_yy_yz_y, g_xx_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_yy_y_x[i] = -2.0 * g_0_yy_yz_x[i] * c_exps[i] + 4.0 * g_xx_yy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_y_y[i] = -2.0 * g_0_yy_yz_y[i] * c_exps[i] + 4.0 * g_xx_yy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_y_z[i] = -2.0 * g_0_yy_yz_z[i] * c_exps[i] + 4.0 * g_xx_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_0_yy_zz_x, g_0_yy_zz_y, g_0_yy_zz_z, g_x_0_z_0_x_yy_z_x, g_x_0_z_0_x_yy_z_y, g_x_0_z_0_x_yy_z_z, g_xx_yy_0_x, g_xx_yy_0_y, g_xx_yy_0_z, g_xx_yy_zz_x, g_xx_yy_zz_y, g_xx_yy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_yy_z_x[i] = g_0_yy_0_x[i] - 2.0 * g_0_yy_zz_x[i] * c_exps[i] - 2.0 * g_xx_yy_0_x[i] * a_exp + 4.0 * g_xx_yy_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_z_y[i] = g_0_yy_0_y[i] - 2.0 * g_0_yy_zz_y[i] * c_exps[i] - 2.0 * g_xx_yy_0_y[i] * a_exp + 4.0 * g_xx_yy_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_z_z[i] = g_0_yy_0_z[i] - 2.0 * g_0_yy_zz_z[i] * c_exps[i] - 2.0 * g_xx_yy_0_z[i] * a_exp + 4.0 * g_xx_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_0_yz_xz_x, g_0_yz_xz_y, g_0_yz_xz_z, g_x_0_z_0_x_yz_x_x, g_x_0_z_0_x_yz_x_y, g_x_0_z_0_x_yz_x_z, g_xx_yz_xz_x, g_xx_yz_xz_y, g_xx_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_yz_x_x[i] = -2.0 * g_0_yz_xz_x[i] * c_exps[i] + 4.0 * g_xx_yz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_x_y[i] = -2.0 * g_0_yz_xz_y[i] * c_exps[i] + 4.0 * g_xx_yz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_x_z[i] = -2.0 * g_0_yz_xz_z[i] * c_exps[i] + 4.0 * g_xx_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_0_yz_yz_x, g_0_yz_yz_y, g_0_yz_yz_z, g_x_0_z_0_x_yz_y_x, g_x_0_z_0_x_yz_y_y, g_x_0_z_0_x_yz_y_z, g_xx_yz_yz_x, g_xx_yz_yz_y, g_xx_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_yz_y_x[i] = -2.0 * g_0_yz_yz_x[i] * c_exps[i] + 4.0 * g_xx_yz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_y_y[i] = -2.0 * g_0_yz_yz_y[i] * c_exps[i] + 4.0 * g_xx_yz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_y_z[i] = -2.0 * g_0_yz_yz_z[i] * c_exps[i] + 4.0 * g_xx_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_0_yz_zz_x, g_0_yz_zz_y, g_0_yz_zz_z, g_x_0_z_0_x_yz_z_x, g_x_0_z_0_x_yz_z_y, g_x_0_z_0_x_yz_z_z, g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z, g_xx_yz_zz_x, g_xx_yz_zz_y, g_xx_yz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_yz_z_x[i] = g_0_yz_0_x[i] - 2.0 * g_0_yz_zz_x[i] * c_exps[i] - 2.0 * g_xx_yz_0_x[i] * a_exp + 4.0 * g_xx_yz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_z_y[i] = g_0_yz_0_y[i] - 2.0 * g_0_yz_zz_y[i] * c_exps[i] - 2.0 * g_xx_yz_0_y[i] * a_exp + 4.0 * g_xx_yz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_z_z[i] = g_0_yz_0_z[i] - 2.0 * g_0_yz_zz_z[i] * c_exps[i] - 2.0 * g_xx_yz_0_z[i] * a_exp + 4.0 * g_xx_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_0_zz_xz_x, g_0_zz_xz_y, g_0_zz_xz_z, g_x_0_z_0_x_zz_x_x, g_x_0_z_0_x_zz_x_y, g_x_0_z_0_x_zz_x_z, g_xx_zz_xz_x, g_xx_zz_xz_y, g_xx_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_zz_x_x[i] = -2.0 * g_0_zz_xz_x[i] * c_exps[i] + 4.0 * g_xx_zz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_x_y[i] = -2.0 * g_0_zz_xz_y[i] * c_exps[i] + 4.0 * g_xx_zz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_x_z[i] = -2.0 * g_0_zz_xz_z[i] * c_exps[i] + 4.0 * g_xx_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_0_zz_yz_x, g_0_zz_yz_y, g_0_zz_yz_z, g_x_0_z_0_x_zz_y_x, g_x_0_z_0_x_zz_y_y, g_x_0_z_0_x_zz_y_z, g_xx_zz_yz_x, g_xx_zz_yz_y, g_xx_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_zz_y_x[i] = -2.0 * g_0_zz_yz_x[i] * c_exps[i] + 4.0 * g_xx_zz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_y_y[i] = -2.0 * g_0_zz_yz_y[i] * c_exps[i] + 4.0 * g_xx_zz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_y_z[i] = -2.0 * g_0_zz_yz_z[i] * c_exps[i] + 4.0 * g_xx_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_0_zz_zz_x, g_0_zz_zz_y, g_0_zz_zz_z, g_x_0_z_0_x_zz_z_x, g_x_0_z_0_x_zz_z_y, g_x_0_z_0_x_zz_z_z, g_xx_zz_0_x, g_xx_zz_0_y, g_xx_zz_0_z, g_xx_zz_zz_x, g_xx_zz_zz_y, g_xx_zz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_zz_z_x[i] = g_0_zz_0_x[i] - 2.0 * g_0_zz_zz_x[i] * c_exps[i] - 2.0 * g_xx_zz_0_x[i] * a_exp + 4.0 * g_xx_zz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_z_y[i] = g_0_zz_0_y[i] - 2.0 * g_0_zz_zz_y[i] * c_exps[i] - 2.0 * g_xx_zz_0_y[i] * a_exp + 4.0 * g_xx_zz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_z_z[i] = g_0_zz_0_z[i] - 2.0 * g_0_zz_zz_z[i] * c_exps[i] - 2.0 * g_xx_zz_0_z[i] * a_exp + 4.0 * g_xx_zz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_x_0_z_0_y_xx_x_x, g_x_0_z_0_y_xx_x_y, g_x_0_z_0_y_xx_x_z, g_xy_xx_xz_x, g_xy_xx_xz_y, g_xy_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xx_x_x[i] = 4.0 * g_xy_xx_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_x_y[i] = 4.0 * g_xy_xx_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_x_z[i] = 4.0 * g_xy_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_x_0_z_0_y_xx_y_x, g_x_0_z_0_y_xx_y_y, g_x_0_z_0_y_xx_y_z, g_xy_xx_yz_x, g_xy_xx_yz_y, g_xy_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xx_y_x[i] = 4.0 * g_xy_xx_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_y_y[i] = 4.0 * g_xy_xx_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_y_z[i] = 4.0 * g_xy_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_x_0_z_0_y_xx_z_x, g_x_0_z_0_y_xx_z_y, g_x_0_z_0_y_xx_z_z, g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xy_xx_zz_x, g_xy_xx_zz_y, g_xy_xx_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xx_z_x[i] = -2.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xy_xx_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_z_y[i] = -2.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xy_xx_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_z_z[i] = -2.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xy_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_x_0_z_0_y_xy_x_x, g_x_0_z_0_y_xy_x_y, g_x_0_z_0_y_xy_x_z, g_xy_xy_xz_x, g_xy_xy_xz_y, g_xy_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xy_x_x[i] = 4.0 * g_xy_xy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_x_y[i] = 4.0 * g_xy_xy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_x_z[i] = 4.0 * g_xy_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_x_0_z_0_y_xy_y_x, g_x_0_z_0_y_xy_y_y, g_x_0_z_0_y_xy_y_z, g_xy_xy_yz_x, g_xy_xy_yz_y, g_xy_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xy_y_x[i] = 4.0 * g_xy_xy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_y_y[i] = 4.0 * g_xy_xy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_y_z[i] = 4.0 * g_xy_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_x_0_z_0_y_xy_z_x, g_x_0_z_0_y_xy_z_y, g_x_0_z_0_y_xy_z_z, g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xy_xy_zz_x, g_xy_xy_zz_y, g_xy_xy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xy_z_x[i] = -2.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xy_xy_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_z_y[i] = -2.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xy_xy_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_z_z[i] = -2.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xy_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_x_0_z_0_y_xz_x_x, g_x_0_z_0_y_xz_x_y, g_x_0_z_0_y_xz_x_z, g_xy_xz_xz_x, g_xy_xz_xz_y, g_xy_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xz_x_x[i] = 4.0 * g_xy_xz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_x_y[i] = 4.0 * g_xy_xz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_x_z[i] = 4.0 * g_xy_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_x_0_z_0_y_xz_y_x, g_x_0_z_0_y_xz_y_y, g_x_0_z_0_y_xz_y_z, g_xy_xz_yz_x, g_xy_xz_yz_y, g_xy_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xz_y_x[i] = 4.0 * g_xy_xz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_y_y[i] = 4.0 * g_xy_xz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_y_z[i] = 4.0 * g_xy_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_x_0_z_0_y_xz_z_x, g_x_0_z_0_y_xz_z_y, g_x_0_z_0_y_xz_z_z, g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xy_xz_zz_x, g_xy_xz_zz_y, g_xy_xz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xz_z_x[i] = -2.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xy_xz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_z_y[i] = -2.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xy_xz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_z_z[i] = -2.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xy_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_x_0_z_0_y_yy_x_x, g_x_0_z_0_y_yy_x_y, g_x_0_z_0_y_yy_x_z, g_xy_yy_xz_x, g_xy_yy_xz_y, g_xy_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_yy_x_x[i] = 4.0 * g_xy_yy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_x_y[i] = 4.0 * g_xy_yy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_x_z[i] = 4.0 * g_xy_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_x_0_z_0_y_yy_y_x, g_x_0_z_0_y_yy_y_y, g_x_0_z_0_y_yy_y_z, g_xy_yy_yz_x, g_xy_yy_yz_y, g_xy_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_yy_y_x[i] = 4.0 * g_xy_yy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_y_y[i] = 4.0 * g_xy_yy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_y_z[i] = 4.0 * g_xy_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_x_0_z_0_y_yy_z_x, g_x_0_z_0_y_yy_z_y, g_x_0_z_0_y_yy_z_z, g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xy_yy_zz_x, g_xy_yy_zz_y, g_xy_yy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_yy_z_x[i] = -2.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xy_yy_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_z_y[i] = -2.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xy_yy_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_z_z[i] = -2.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xy_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_x_0_z_0_y_yz_x_x, g_x_0_z_0_y_yz_x_y, g_x_0_z_0_y_yz_x_z, g_xy_yz_xz_x, g_xy_yz_xz_y, g_xy_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_yz_x_x[i] = 4.0 * g_xy_yz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_x_y[i] = 4.0 * g_xy_yz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_x_z[i] = 4.0 * g_xy_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_x_0_z_0_y_yz_y_x, g_x_0_z_0_y_yz_y_y, g_x_0_z_0_y_yz_y_z, g_xy_yz_yz_x, g_xy_yz_yz_y, g_xy_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_yz_y_x[i] = 4.0 * g_xy_yz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_y_y[i] = 4.0 * g_xy_yz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_y_z[i] = 4.0 * g_xy_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_x_0_z_0_y_yz_z_x, g_x_0_z_0_y_yz_z_y, g_x_0_z_0_y_yz_z_z, g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xy_yz_zz_x, g_xy_yz_zz_y, g_xy_yz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_yz_z_x[i] = -2.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xy_yz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_z_y[i] = -2.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xy_yz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_z_z[i] = -2.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xy_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_x_0_z_0_y_zz_x_x, g_x_0_z_0_y_zz_x_y, g_x_0_z_0_y_zz_x_z, g_xy_zz_xz_x, g_xy_zz_xz_y, g_xy_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_zz_x_x[i] = 4.0 * g_xy_zz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_x_y[i] = 4.0 * g_xy_zz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_x_z[i] = 4.0 * g_xy_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_x_0_z_0_y_zz_y_x, g_x_0_z_0_y_zz_y_y, g_x_0_z_0_y_zz_y_z, g_xy_zz_yz_x, g_xy_zz_yz_y, g_xy_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_zz_y_x[i] = 4.0 * g_xy_zz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_y_y[i] = 4.0 * g_xy_zz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_y_z[i] = 4.0 * g_xy_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_x_0_z_0_y_zz_z_x, g_x_0_z_0_y_zz_z_y, g_x_0_z_0_y_zz_z_z, g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xy_zz_zz_x, g_xy_zz_zz_y, g_xy_zz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_zz_z_x[i] = -2.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xy_zz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_z_y[i] = -2.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xy_zz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_z_z[i] = -2.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xy_zz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_x_0_z_0_z_xx_x_x, g_x_0_z_0_z_xx_x_y, g_x_0_z_0_z_xx_x_z, g_xz_xx_xz_x, g_xz_xx_xz_y, g_xz_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xx_x_x[i] = 4.0 * g_xz_xx_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_x_y[i] = 4.0 * g_xz_xx_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_x_z[i] = 4.0 * g_xz_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_x_0_z_0_z_xx_y_x, g_x_0_z_0_z_xx_y_y, g_x_0_z_0_z_xx_y_z, g_xz_xx_yz_x, g_xz_xx_yz_y, g_xz_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xx_y_x[i] = 4.0 * g_xz_xx_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_y_y[i] = 4.0 * g_xz_xx_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_y_z[i] = 4.0 * g_xz_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_x_0_z_0_z_xx_z_x, g_x_0_z_0_z_xx_z_y, g_x_0_z_0_z_xx_z_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_xz_xx_zz_x, g_xz_xx_zz_y, g_xz_xx_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xx_z_x[i] = -2.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xz_xx_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_z_y[i] = -2.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xz_xx_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_z_z[i] = -2.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xz_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_x_0_z_0_z_xy_x_x, g_x_0_z_0_z_xy_x_y, g_x_0_z_0_z_xy_x_z, g_xz_xy_xz_x, g_xz_xy_xz_y, g_xz_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xy_x_x[i] = 4.0 * g_xz_xy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_x_y[i] = 4.0 * g_xz_xy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_x_z[i] = 4.0 * g_xz_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_x_0_z_0_z_xy_y_x, g_x_0_z_0_z_xy_y_y, g_x_0_z_0_z_xy_y_z, g_xz_xy_yz_x, g_xz_xy_yz_y, g_xz_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xy_y_x[i] = 4.0 * g_xz_xy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_y_y[i] = 4.0 * g_xz_xy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_y_z[i] = 4.0 * g_xz_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_x_0_z_0_z_xy_z_x, g_x_0_z_0_z_xy_z_y, g_x_0_z_0_z_xy_z_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_xz_xy_zz_x, g_xz_xy_zz_y, g_xz_xy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xy_z_x[i] = -2.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xz_xy_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_z_y[i] = -2.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xz_xy_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_z_z[i] = -2.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xz_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_x_0_z_0_z_xz_x_x, g_x_0_z_0_z_xz_x_y, g_x_0_z_0_z_xz_x_z, g_xz_xz_xz_x, g_xz_xz_xz_y, g_xz_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xz_x_x[i] = 4.0 * g_xz_xz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_x_y[i] = 4.0 * g_xz_xz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_x_z[i] = 4.0 * g_xz_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_x_0_z_0_z_xz_y_x, g_x_0_z_0_z_xz_y_y, g_x_0_z_0_z_xz_y_z, g_xz_xz_yz_x, g_xz_xz_yz_y, g_xz_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xz_y_x[i] = 4.0 * g_xz_xz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_y_y[i] = 4.0 * g_xz_xz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_y_z[i] = 4.0 * g_xz_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_x_0_z_0_z_xz_z_x, g_x_0_z_0_z_xz_z_y, g_x_0_z_0_z_xz_z_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_xz_xz_zz_x, g_xz_xz_zz_y, g_xz_xz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xz_z_x[i] = -2.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xz_xz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_z_y[i] = -2.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xz_xz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_z_z[i] = -2.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xz_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_x_0_z_0_z_yy_x_x, g_x_0_z_0_z_yy_x_y, g_x_0_z_0_z_yy_x_z, g_xz_yy_xz_x, g_xz_yy_xz_y, g_xz_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_yy_x_x[i] = 4.0 * g_xz_yy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_x_y[i] = 4.0 * g_xz_yy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_x_z[i] = 4.0 * g_xz_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_x_0_z_0_z_yy_y_x, g_x_0_z_0_z_yy_y_y, g_x_0_z_0_z_yy_y_z, g_xz_yy_yz_x, g_xz_yy_yz_y, g_xz_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_yy_y_x[i] = 4.0 * g_xz_yy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_y_y[i] = 4.0 * g_xz_yy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_y_z[i] = 4.0 * g_xz_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_x_0_z_0_z_yy_z_x, g_x_0_z_0_z_yy_z_y, g_x_0_z_0_z_yy_z_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_xz_yy_zz_x, g_xz_yy_zz_y, g_xz_yy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_yy_z_x[i] = -2.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xz_yy_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_z_y[i] = -2.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xz_yy_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_z_z[i] = -2.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xz_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_x_0_z_0_z_yz_x_x, g_x_0_z_0_z_yz_x_y, g_x_0_z_0_z_yz_x_z, g_xz_yz_xz_x, g_xz_yz_xz_y, g_xz_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_yz_x_x[i] = 4.0 * g_xz_yz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_x_y[i] = 4.0 * g_xz_yz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_x_z[i] = 4.0 * g_xz_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_x_0_z_0_z_yz_y_x, g_x_0_z_0_z_yz_y_y, g_x_0_z_0_z_yz_y_z, g_xz_yz_yz_x, g_xz_yz_yz_y, g_xz_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_yz_y_x[i] = 4.0 * g_xz_yz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_y_y[i] = 4.0 * g_xz_yz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_y_z[i] = 4.0 * g_xz_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_x_0_z_0_z_yz_z_x, g_x_0_z_0_z_yz_z_y, g_x_0_z_0_z_yz_z_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_xz_yz_zz_x, g_xz_yz_zz_y, g_xz_yz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_yz_z_x[i] = -2.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xz_yz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_z_y[i] = -2.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xz_yz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_z_z[i] = -2.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xz_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_x_0_z_0_z_zz_x_x, g_x_0_z_0_z_zz_x_y, g_x_0_z_0_z_zz_x_z, g_xz_zz_xz_x, g_xz_zz_xz_y, g_xz_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_zz_x_x[i] = 4.0 * g_xz_zz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_x_y[i] = 4.0 * g_xz_zz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_x_z[i] = 4.0 * g_xz_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_x_0_z_0_z_zz_y_x, g_x_0_z_0_z_zz_y_y, g_x_0_z_0_z_zz_y_z, g_xz_zz_yz_x, g_xz_zz_yz_y, g_xz_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_zz_y_x[i] = 4.0 * g_xz_zz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_y_y[i] = 4.0 * g_xz_zz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_y_z[i] = 4.0 * g_xz_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_x_0_z_0_z_zz_z_x, g_x_0_z_0_z_zz_z_y, g_x_0_z_0_z_zz_z_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_xz_zz_zz_x, g_xz_zz_zz_y, g_xz_zz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_zz_z_x[i] = -2.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xz_zz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_z_y[i] = -2.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xz_zz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_z_z[i] = -2.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xz_zz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (486-489)

    #pragma omp simd aligned(g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xy_xx_xx_x, g_xy_xx_xx_y, g_xy_xx_xx_z, g_y_0_x_0_x_xx_x_x, g_y_0_x_0_x_xx_x_y, g_y_0_x_0_x_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xx_x_x[i] = -2.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xy_xx_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_x_y[i] = -2.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xy_xx_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_x_z[i] = -2.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xy_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (489-492)

    #pragma omp simd aligned(g_xy_xx_xy_x, g_xy_xx_xy_y, g_xy_xx_xy_z, g_y_0_x_0_x_xx_y_x, g_y_0_x_0_x_xx_y_y, g_y_0_x_0_x_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xx_y_x[i] = 4.0 * g_xy_xx_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_y_y[i] = 4.0 * g_xy_xx_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_y_z[i] = 4.0 * g_xy_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (492-495)

    #pragma omp simd aligned(g_xy_xx_xz_x, g_xy_xx_xz_y, g_xy_xx_xz_z, g_y_0_x_0_x_xx_z_x, g_y_0_x_0_x_xx_z_y, g_y_0_x_0_x_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xx_z_x[i] = 4.0 * g_xy_xx_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_z_y[i] = 4.0 * g_xy_xx_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_z_z[i] = 4.0 * g_xy_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (495-498)

    #pragma omp simd aligned(g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xy_xy_xx_x, g_xy_xy_xx_y, g_xy_xy_xx_z, g_y_0_x_0_x_xy_x_x, g_y_0_x_0_x_xy_x_y, g_y_0_x_0_x_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xy_x_x[i] = -2.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xy_xy_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_x_y[i] = -2.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xy_xy_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_x_z[i] = -2.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xy_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (498-501)

    #pragma omp simd aligned(g_xy_xy_xy_x, g_xy_xy_xy_y, g_xy_xy_xy_z, g_y_0_x_0_x_xy_y_x, g_y_0_x_0_x_xy_y_y, g_y_0_x_0_x_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xy_y_x[i] = 4.0 * g_xy_xy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_y_y[i] = 4.0 * g_xy_xy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_y_z[i] = 4.0 * g_xy_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (501-504)

    #pragma omp simd aligned(g_xy_xy_xz_x, g_xy_xy_xz_y, g_xy_xy_xz_z, g_y_0_x_0_x_xy_z_x, g_y_0_x_0_x_xy_z_y, g_y_0_x_0_x_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xy_z_x[i] = 4.0 * g_xy_xy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_z_y[i] = 4.0 * g_xy_xy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_z_z[i] = 4.0 * g_xy_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (504-507)

    #pragma omp simd aligned(g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xy_xz_xx_x, g_xy_xz_xx_y, g_xy_xz_xx_z, g_y_0_x_0_x_xz_x_x, g_y_0_x_0_x_xz_x_y, g_y_0_x_0_x_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xz_x_x[i] = -2.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xy_xz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_x_y[i] = -2.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xy_xz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_x_z[i] = -2.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xy_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (507-510)

    #pragma omp simd aligned(g_xy_xz_xy_x, g_xy_xz_xy_y, g_xy_xz_xy_z, g_y_0_x_0_x_xz_y_x, g_y_0_x_0_x_xz_y_y, g_y_0_x_0_x_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xz_y_x[i] = 4.0 * g_xy_xz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_y_y[i] = 4.0 * g_xy_xz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_y_z[i] = 4.0 * g_xy_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (510-513)

    #pragma omp simd aligned(g_xy_xz_xz_x, g_xy_xz_xz_y, g_xy_xz_xz_z, g_y_0_x_0_x_xz_z_x, g_y_0_x_0_x_xz_z_y, g_y_0_x_0_x_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xz_z_x[i] = 4.0 * g_xy_xz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_z_y[i] = 4.0 * g_xy_xz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_z_z[i] = 4.0 * g_xy_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (513-516)

    #pragma omp simd aligned(g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xy_yy_xx_x, g_xy_yy_xx_y, g_xy_yy_xx_z, g_y_0_x_0_x_yy_x_x, g_y_0_x_0_x_yy_x_y, g_y_0_x_0_x_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_yy_x_x[i] = -2.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xy_yy_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_x_y[i] = -2.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xy_yy_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_x_z[i] = -2.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xy_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (516-519)

    #pragma omp simd aligned(g_xy_yy_xy_x, g_xy_yy_xy_y, g_xy_yy_xy_z, g_y_0_x_0_x_yy_y_x, g_y_0_x_0_x_yy_y_y, g_y_0_x_0_x_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_yy_y_x[i] = 4.0 * g_xy_yy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_y_y[i] = 4.0 * g_xy_yy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_y_z[i] = 4.0 * g_xy_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (519-522)

    #pragma omp simd aligned(g_xy_yy_xz_x, g_xy_yy_xz_y, g_xy_yy_xz_z, g_y_0_x_0_x_yy_z_x, g_y_0_x_0_x_yy_z_y, g_y_0_x_0_x_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_yy_z_x[i] = 4.0 * g_xy_yy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_z_y[i] = 4.0 * g_xy_yy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_z_z[i] = 4.0 * g_xy_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (522-525)

    #pragma omp simd aligned(g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xy_yz_xx_x, g_xy_yz_xx_y, g_xy_yz_xx_z, g_y_0_x_0_x_yz_x_x, g_y_0_x_0_x_yz_x_y, g_y_0_x_0_x_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_yz_x_x[i] = -2.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xy_yz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_x_y[i] = -2.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xy_yz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_x_z[i] = -2.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xy_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (525-528)

    #pragma omp simd aligned(g_xy_yz_xy_x, g_xy_yz_xy_y, g_xy_yz_xy_z, g_y_0_x_0_x_yz_y_x, g_y_0_x_0_x_yz_y_y, g_y_0_x_0_x_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_yz_y_x[i] = 4.0 * g_xy_yz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_y_y[i] = 4.0 * g_xy_yz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_y_z[i] = 4.0 * g_xy_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (528-531)

    #pragma omp simd aligned(g_xy_yz_xz_x, g_xy_yz_xz_y, g_xy_yz_xz_z, g_y_0_x_0_x_yz_z_x, g_y_0_x_0_x_yz_z_y, g_y_0_x_0_x_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_yz_z_x[i] = 4.0 * g_xy_yz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_z_y[i] = 4.0 * g_xy_yz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_z_z[i] = 4.0 * g_xy_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (531-534)

    #pragma omp simd aligned(g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xy_zz_xx_x, g_xy_zz_xx_y, g_xy_zz_xx_z, g_y_0_x_0_x_zz_x_x, g_y_0_x_0_x_zz_x_y, g_y_0_x_0_x_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_zz_x_x[i] = -2.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xy_zz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_x_y[i] = -2.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xy_zz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_x_z[i] = -2.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xy_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (534-537)

    #pragma omp simd aligned(g_xy_zz_xy_x, g_xy_zz_xy_y, g_xy_zz_xy_z, g_y_0_x_0_x_zz_y_x, g_y_0_x_0_x_zz_y_y, g_y_0_x_0_x_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_zz_y_x[i] = 4.0 * g_xy_zz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_y_y[i] = 4.0 * g_xy_zz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_y_z[i] = 4.0 * g_xy_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (537-540)

    #pragma omp simd aligned(g_xy_zz_xz_x, g_xy_zz_xz_y, g_xy_zz_xz_z, g_y_0_x_0_x_zz_z_x, g_y_0_x_0_x_zz_z_y, g_y_0_x_0_x_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_zz_z_x[i] = 4.0 * g_xy_zz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_z_y[i] = 4.0 * g_xy_zz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_z_z[i] = 4.0 * g_xy_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (540-543)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_0_xx_xx_x, g_0_xx_xx_y, g_0_xx_xx_z, g_y_0_x_0_y_xx_x_x, g_y_0_x_0_y_xx_x_y, g_y_0_x_0_y_xx_x_z, g_yy_xx_0_x, g_yy_xx_0_y, g_yy_xx_0_z, g_yy_xx_xx_x, g_yy_xx_xx_y, g_yy_xx_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xx_x_x[i] = g_0_xx_0_x[i] - 2.0 * g_0_xx_xx_x[i] * c_exps[i] - 2.0 * g_yy_xx_0_x[i] * a_exp + 4.0 * g_yy_xx_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_x_y[i] = g_0_xx_0_y[i] - 2.0 * g_0_xx_xx_y[i] * c_exps[i] - 2.0 * g_yy_xx_0_y[i] * a_exp + 4.0 * g_yy_xx_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_x_z[i] = g_0_xx_0_z[i] - 2.0 * g_0_xx_xx_z[i] * c_exps[i] - 2.0 * g_yy_xx_0_z[i] * a_exp + 4.0 * g_yy_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (543-546)

    #pragma omp simd aligned(g_0_xx_xy_x, g_0_xx_xy_y, g_0_xx_xy_z, g_y_0_x_0_y_xx_y_x, g_y_0_x_0_y_xx_y_y, g_y_0_x_0_y_xx_y_z, g_yy_xx_xy_x, g_yy_xx_xy_y, g_yy_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xx_y_x[i] = -2.0 * g_0_xx_xy_x[i] * c_exps[i] + 4.0 * g_yy_xx_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_y_y[i] = -2.0 * g_0_xx_xy_y[i] * c_exps[i] + 4.0 * g_yy_xx_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_y_z[i] = -2.0 * g_0_xx_xy_z[i] * c_exps[i] + 4.0 * g_yy_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (546-549)

    #pragma omp simd aligned(g_0_xx_xz_x, g_0_xx_xz_y, g_0_xx_xz_z, g_y_0_x_0_y_xx_z_x, g_y_0_x_0_y_xx_z_y, g_y_0_x_0_y_xx_z_z, g_yy_xx_xz_x, g_yy_xx_xz_y, g_yy_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xx_z_x[i] = -2.0 * g_0_xx_xz_x[i] * c_exps[i] + 4.0 * g_yy_xx_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_z_y[i] = -2.0 * g_0_xx_xz_y[i] * c_exps[i] + 4.0 * g_yy_xx_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_z_z[i] = -2.0 * g_0_xx_xz_z[i] * c_exps[i] + 4.0 * g_yy_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (549-552)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_0_xy_xx_x, g_0_xy_xx_y, g_0_xy_xx_z, g_y_0_x_0_y_xy_x_x, g_y_0_x_0_y_xy_x_y, g_y_0_x_0_y_xy_x_z, g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z, g_yy_xy_xx_x, g_yy_xy_xx_y, g_yy_xy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xy_x_x[i] = g_0_xy_0_x[i] - 2.0 * g_0_xy_xx_x[i] * c_exps[i] - 2.0 * g_yy_xy_0_x[i] * a_exp + 4.0 * g_yy_xy_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_x_y[i] = g_0_xy_0_y[i] - 2.0 * g_0_xy_xx_y[i] * c_exps[i] - 2.0 * g_yy_xy_0_y[i] * a_exp + 4.0 * g_yy_xy_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_x_z[i] = g_0_xy_0_z[i] - 2.0 * g_0_xy_xx_z[i] * c_exps[i] - 2.0 * g_yy_xy_0_z[i] * a_exp + 4.0 * g_yy_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (552-555)

    #pragma omp simd aligned(g_0_xy_xy_x, g_0_xy_xy_y, g_0_xy_xy_z, g_y_0_x_0_y_xy_y_x, g_y_0_x_0_y_xy_y_y, g_y_0_x_0_y_xy_y_z, g_yy_xy_xy_x, g_yy_xy_xy_y, g_yy_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xy_y_x[i] = -2.0 * g_0_xy_xy_x[i] * c_exps[i] + 4.0 * g_yy_xy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_y_y[i] = -2.0 * g_0_xy_xy_y[i] * c_exps[i] + 4.0 * g_yy_xy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_y_z[i] = -2.0 * g_0_xy_xy_z[i] * c_exps[i] + 4.0 * g_yy_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (555-558)

    #pragma omp simd aligned(g_0_xy_xz_x, g_0_xy_xz_y, g_0_xy_xz_z, g_y_0_x_0_y_xy_z_x, g_y_0_x_0_y_xy_z_y, g_y_0_x_0_y_xy_z_z, g_yy_xy_xz_x, g_yy_xy_xz_y, g_yy_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xy_z_x[i] = -2.0 * g_0_xy_xz_x[i] * c_exps[i] + 4.0 * g_yy_xy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_z_y[i] = -2.0 * g_0_xy_xz_y[i] * c_exps[i] + 4.0 * g_yy_xy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_z_z[i] = -2.0 * g_0_xy_xz_z[i] * c_exps[i] + 4.0 * g_yy_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (558-561)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_0_xz_xx_x, g_0_xz_xx_y, g_0_xz_xx_z, g_y_0_x_0_y_xz_x_x, g_y_0_x_0_y_xz_x_y, g_y_0_x_0_y_xz_x_z, g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z, g_yy_xz_xx_x, g_yy_xz_xx_y, g_yy_xz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xz_x_x[i] = g_0_xz_0_x[i] - 2.0 * g_0_xz_xx_x[i] * c_exps[i] - 2.0 * g_yy_xz_0_x[i] * a_exp + 4.0 * g_yy_xz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_x_y[i] = g_0_xz_0_y[i] - 2.0 * g_0_xz_xx_y[i] * c_exps[i] - 2.0 * g_yy_xz_0_y[i] * a_exp + 4.0 * g_yy_xz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_x_z[i] = g_0_xz_0_z[i] - 2.0 * g_0_xz_xx_z[i] * c_exps[i] - 2.0 * g_yy_xz_0_z[i] * a_exp + 4.0 * g_yy_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (561-564)

    #pragma omp simd aligned(g_0_xz_xy_x, g_0_xz_xy_y, g_0_xz_xy_z, g_y_0_x_0_y_xz_y_x, g_y_0_x_0_y_xz_y_y, g_y_0_x_0_y_xz_y_z, g_yy_xz_xy_x, g_yy_xz_xy_y, g_yy_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xz_y_x[i] = -2.0 * g_0_xz_xy_x[i] * c_exps[i] + 4.0 * g_yy_xz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_y_y[i] = -2.0 * g_0_xz_xy_y[i] * c_exps[i] + 4.0 * g_yy_xz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_y_z[i] = -2.0 * g_0_xz_xy_z[i] * c_exps[i] + 4.0 * g_yy_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (564-567)

    #pragma omp simd aligned(g_0_xz_xz_x, g_0_xz_xz_y, g_0_xz_xz_z, g_y_0_x_0_y_xz_z_x, g_y_0_x_0_y_xz_z_y, g_y_0_x_0_y_xz_z_z, g_yy_xz_xz_x, g_yy_xz_xz_y, g_yy_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xz_z_x[i] = -2.0 * g_0_xz_xz_x[i] * c_exps[i] + 4.0 * g_yy_xz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_z_y[i] = -2.0 * g_0_xz_xz_y[i] * c_exps[i] + 4.0 * g_yy_xz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_z_z[i] = -2.0 * g_0_xz_xz_z[i] * c_exps[i] + 4.0 * g_yy_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (567-570)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_0_yy_xx_x, g_0_yy_xx_y, g_0_yy_xx_z, g_y_0_x_0_y_yy_x_x, g_y_0_x_0_y_yy_x_y, g_y_0_x_0_y_yy_x_z, g_yy_yy_0_x, g_yy_yy_0_y, g_yy_yy_0_z, g_yy_yy_xx_x, g_yy_yy_xx_y, g_yy_yy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_yy_x_x[i] = g_0_yy_0_x[i] - 2.0 * g_0_yy_xx_x[i] * c_exps[i] - 2.0 * g_yy_yy_0_x[i] * a_exp + 4.0 * g_yy_yy_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_x_y[i] = g_0_yy_0_y[i] - 2.0 * g_0_yy_xx_y[i] * c_exps[i] - 2.0 * g_yy_yy_0_y[i] * a_exp + 4.0 * g_yy_yy_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_x_z[i] = g_0_yy_0_z[i] - 2.0 * g_0_yy_xx_z[i] * c_exps[i] - 2.0 * g_yy_yy_0_z[i] * a_exp + 4.0 * g_yy_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (570-573)

    #pragma omp simd aligned(g_0_yy_xy_x, g_0_yy_xy_y, g_0_yy_xy_z, g_y_0_x_0_y_yy_y_x, g_y_0_x_0_y_yy_y_y, g_y_0_x_0_y_yy_y_z, g_yy_yy_xy_x, g_yy_yy_xy_y, g_yy_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_yy_y_x[i] = -2.0 * g_0_yy_xy_x[i] * c_exps[i] + 4.0 * g_yy_yy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_y_y[i] = -2.0 * g_0_yy_xy_y[i] * c_exps[i] + 4.0 * g_yy_yy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_y_z[i] = -2.0 * g_0_yy_xy_z[i] * c_exps[i] + 4.0 * g_yy_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (573-576)

    #pragma omp simd aligned(g_0_yy_xz_x, g_0_yy_xz_y, g_0_yy_xz_z, g_y_0_x_0_y_yy_z_x, g_y_0_x_0_y_yy_z_y, g_y_0_x_0_y_yy_z_z, g_yy_yy_xz_x, g_yy_yy_xz_y, g_yy_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_yy_z_x[i] = -2.0 * g_0_yy_xz_x[i] * c_exps[i] + 4.0 * g_yy_yy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_z_y[i] = -2.0 * g_0_yy_xz_y[i] * c_exps[i] + 4.0 * g_yy_yy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_z_z[i] = -2.0 * g_0_yy_xz_z[i] * c_exps[i] + 4.0 * g_yy_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (576-579)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_0_yz_xx_x, g_0_yz_xx_y, g_0_yz_xx_z, g_y_0_x_0_y_yz_x_x, g_y_0_x_0_y_yz_x_y, g_y_0_x_0_y_yz_x_z, g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z, g_yy_yz_xx_x, g_yy_yz_xx_y, g_yy_yz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_yz_x_x[i] = g_0_yz_0_x[i] - 2.0 * g_0_yz_xx_x[i] * c_exps[i] - 2.0 * g_yy_yz_0_x[i] * a_exp + 4.0 * g_yy_yz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_x_y[i] = g_0_yz_0_y[i] - 2.0 * g_0_yz_xx_y[i] * c_exps[i] - 2.0 * g_yy_yz_0_y[i] * a_exp + 4.0 * g_yy_yz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_x_z[i] = g_0_yz_0_z[i] - 2.0 * g_0_yz_xx_z[i] * c_exps[i] - 2.0 * g_yy_yz_0_z[i] * a_exp + 4.0 * g_yy_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (579-582)

    #pragma omp simd aligned(g_0_yz_xy_x, g_0_yz_xy_y, g_0_yz_xy_z, g_y_0_x_0_y_yz_y_x, g_y_0_x_0_y_yz_y_y, g_y_0_x_0_y_yz_y_z, g_yy_yz_xy_x, g_yy_yz_xy_y, g_yy_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_yz_y_x[i] = -2.0 * g_0_yz_xy_x[i] * c_exps[i] + 4.0 * g_yy_yz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_y_y[i] = -2.0 * g_0_yz_xy_y[i] * c_exps[i] + 4.0 * g_yy_yz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_y_z[i] = -2.0 * g_0_yz_xy_z[i] * c_exps[i] + 4.0 * g_yy_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (582-585)

    #pragma omp simd aligned(g_0_yz_xz_x, g_0_yz_xz_y, g_0_yz_xz_z, g_y_0_x_0_y_yz_z_x, g_y_0_x_0_y_yz_z_y, g_y_0_x_0_y_yz_z_z, g_yy_yz_xz_x, g_yy_yz_xz_y, g_yy_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_yz_z_x[i] = -2.0 * g_0_yz_xz_x[i] * c_exps[i] + 4.0 * g_yy_yz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_z_y[i] = -2.0 * g_0_yz_xz_y[i] * c_exps[i] + 4.0 * g_yy_yz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_z_z[i] = -2.0 * g_0_yz_xz_z[i] * c_exps[i] + 4.0 * g_yy_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (585-588)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_0_zz_xx_x, g_0_zz_xx_y, g_0_zz_xx_z, g_y_0_x_0_y_zz_x_x, g_y_0_x_0_y_zz_x_y, g_y_0_x_0_y_zz_x_z, g_yy_zz_0_x, g_yy_zz_0_y, g_yy_zz_0_z, g_yy_zz_xx_x, g_yy_zz_xx_y, g_yy_zz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_zz_x_x[i] = g_0_zz_0_x[i] - 2.0 * g_0_zz_xx_x[i] * c_exps[i] - 2.0 * g_yy_zz_0_x[i] * a_exp + 4.0 * g_yy_zz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_x_y[i] = g_0_zz_0_y[i] - 2.0 * g_0_zz_xx_y[i] * c_exps[i] - 2.0 * g_yy_zz_0_y[i] * a_exp + 4.0 * g_yy_zz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_x_z[i] = g_0_zz_0_z[i] - 2.0 * g_0_zz_xx_z[i] * c_exps[i] - 2.0 * g_yy_zz_0_z[i] * a_exp + 4.0 * g_yy_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (588-591)

    #pragma omp simd aligned(g_0_zz_xy_x, g_0_zz_xy_y, g_0_zz_xy_z, g_y_0_x_0_y_zz_y_x, g_y_0_x_0_y_zz_y_y, g_y_0_x_0_y_zz_y_z, g_yy_zz_xy_x, g_yy_zz_xy_y, g_yy_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_zz_y_x[i] = -2.0 * g_0_zz_xy_x[i] * c_exps[i] + 4.0 * g_yy_zz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_y_y[i] = -2.0 * g_0_zz_xy_y[i] * c_exps[i] + 4.0 * g_yy_zz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_y_z[i] = -2.0 * g_0_zz_xy_z[i] * c_exps[i] + 4.0 * g_yy_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (591-594)

    #pragma omp simd aligned(g_0_zz_xz_x, g_0_zz_xz_y, g_0_zz_xz_z, g_y_0_x_0_y_zz_z_x, g_y_0_x_0_y_zz_z_y, g_y_0_x_0_y_zz_z_z, g_yy_zz_xz_x, g_yy_zz_xz_y, g_yy_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_zz_z_x[i] = -2.0 * g_0_zz_xz_x[i] * c_exps[i] + 4.0 * g_yy_zz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_z_y[i] = -2.0 * g_0_zz_xz_y[i] * c_exps[i] + 4.0 * g_yy_zz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_z_z[i] = -2.0 * g_0_zz_xz_z[i] * c_exps[i] + 4.0 * g_yy_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (594-597)

    #pragma omp simd aligned(g_y_0_x_0_z_xx_x_x, g_y_0_x_0_z_xx_x_y, g_y_0_x_0_z_xx_x_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z, g_yz_xx_xx_x, g_yz_xx_xx_y, g_yz_xx_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xx_x_x[i] = -2.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_yz_xx_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_x_y[i] = -2.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_yz_xx_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_x_z[i] = -2.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_yz_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (597-600)

    #pragma omp simd aligned(g_y_0_x_0_z_xx_y_x, g_y_0_x_0_z_xx_y_y, g_y_0_x_0_z_xx_y_z, g_yz_xx_xy_x, g_yz_xx_xy_y, g_yz_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xx_y_x[i] = 4.0 * g_yz_xx_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_y_y[i] = 4.0 * g_yz_xx_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_y_z[i] = 4.0 * g_yz_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (600-603)

    #pragma omp simd aligned(g_y_0_x_0_z_xx_z_x, g_y_0_x_0_z_xx_z_y, g_y_0_x_0_z_xx_z_z, g_yz_xx_xz_x, g_yz_xx_xz_y, g_yz_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xx_z_x[i] = 4.0 * g_yz_xx_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_z_y[i] = 4.0 * g_yz_xx_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_z_z[i] = 4.0 * g_yz_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (603-606)

    #pragma omp simd aligned(g_y_0_x_0_z_xy_x_x, g_y_0_x_0_z_xy_x_y, g_y_0_x_0_z_xy_x_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_yz_xy_xx_x, g_yz_xy_xx_y, g_yz_xy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xy_x_x[i] = -2.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_yz_xy_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_x_y[i] = -2.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_yz_xy_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_x_z[i] = -2.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_yz_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (606-609)

    #pragma omp simd aligned(g_y_0_x_0_z_xy_y_x, g_y_0_x_0_z_xy_y_y, g_y_0_x_0_z_xy_y_z, g_yz_xy_xy_x, g_yz_xy_xy_y, g_yz_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xy_y_x[i] = 4.0 * g_yz_xy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_y_y[i] = 4.0 * g_yz_xy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_y_z[i] = 4.0 * g_yz_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (609-612)

    #pragma omp simd aligned(g_y_0_x_0_z_xy_z_x, g_y_0_x_0_z_xy_z_y, g_y_0_x_0_z_xy_z_z, g_yz_xy_xz_x, g_yz_xy_xz_y, g_yz_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xy_z_x[i] = 4.0 * g_yz_xy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_z_y[i] = 4.0 * g_yz_xy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_z_z[i] = 4.0 * g_yz_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (612-615)

    #pragma omp simd aligned(g_y_0_x_0_z_xz_x_x, g_y_0_x_0_z_xz_x_y, g_y_0_x_0_z_xz_x_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_yz_xz_xx_x, g_yz_xz_xx_y, g_yz_xz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xz_x_x[i] = -2.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_yz_xz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_x_y[i] = -2.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_yz_xz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_x_z[i] = -2.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_yz_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (615-618)

    #pragma omp simd aligned(g_y_0_x_0_z_xz_y_x, g_y_0_x_0_z_xz_y_y, g_y_0_x_0_z_xz_y_z, g_yz_xz_xy_x, g_yz_xz_xy_y, g_yz_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xz_y_x[i] = 4.0 * g_yz_xz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_y_y[i] = 4.0 * g_yz_xz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_y_z[i] = 4.0 * g_yz_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (618-621)

    #pragma omp simd aligned(g_y_0_x_0_z_xz_z_x, g_y_0_x_0_z_xz_z_y, g_y_0_x_0_z_xz_z_z, g_yz_xz_xz_x, g_yz_xz_xz_y, g_yz_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xz_z_x[i] = 4.0 * g_yz_xz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_z_y[i] = 4.0 * g_yz_xz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_z_z[i] = 4.0 * g_yz_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (621-624)

    #pragma omp simd aligned(g_y_0_x_0_z_yy_x_x, g_y_0_x_0_z_yy_x_y, g_y_0_x_0_z_yy_x_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z, g_yz_yy_xx_x, g_yz_yy_xx_y, g_yz_yy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_yy_x_x[i] = -2.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_yz_yy_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_x_y[i] = -2.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_yz_yy_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_x_z[i] = -2.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_yz_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (624-627)

    #pragma omp simd aligned(g_y_0_x_0_z_yy_y_x, g_y_0_x_0_z_yy_y_y, g_y_0_x_0_z_yy_y_z, g_yz_yy_xy_x, g_yz_yy_xy_y, g_yz_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_yy_y_x[i] = 4.0 * g_yz_yy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_y_y[i] = 4.0 * g_yz_yy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_y_z[i] = 4.0 * g_yz_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (627-630)

    #pragma omp simd aligned(g_y_0_x_0_z_yy_z_x, g_y_0_x_0_z_yy_z_y, g_y_0_x_0_z_yy_z_z, g_yz_yy_xz_x, g_yz_yy_xz_y, g_yz_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_yy_z_x[i] = 4.0 * g_yz_yy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_z_y[i] = 4.0 * g_yz_yy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_z_z[i] = 4.0 * g_yz_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (630-633)

    #pragma omp simd aligned(g_y_0_x_0_z_yz_x_x, g_y_0_x_0_z_yz_x_y, g_y_0_x_0_z_yz_x_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_yz_yz_xx_x, g_yz_yz_xx_y, g_yz_yz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_yz_x_x[i] = -2.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_yz_yz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_x_y[i] = -2.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_yz_yz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_x_z[i] = -2.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_yz_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (633-636)

    #pragma omp simd aligned(g_y_0_x_0_z_yz_y_x, g_y_0_x_0_z_yz_y_y, g_y_0_x_0_z_yz_y_z, g_yz_yz_xy_x, g_yz_yz_xy_y, g_yz_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_yz_y_x[i] = 4.0 * g_yz_yz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_y_y[i] = 4.0 * g_yz_yz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_y_z[i] = 4.0 * g_yz_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (636-639)

    #pragma omp simd aligned(g_y_0_x_0_z_yz_z_x, g_y_0_x_0_z_yz_z_y, g_y_0_x_0_z_yz_z_z, g_yz_yz_xz_x, g_yz_yz_xz_y, g_yz_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_yz_z_x[i] = 4.0 * g_yz_yz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_z_y[i] = 4.0 * g_yz_yz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_z_z[i] = 4.0 * g_yz_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (639-642)

    #pragma omp simd aligned(g_y_0_x_0_z_zz_x_x, g_y_0_x_0_z_zz_x_y, g_y_0_x_0_z_zz_x_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z, g_yz_zz_xx_x, g_yz_zz_xx_y, g_yz_zz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_zz_x_x[i] = -2.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_yz_zz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_x_y[i] = -2.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_yz_zz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_x_z[i] = -2.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_yz_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (642-645)

    #pragma omp simd aligned(g_y_0_x_0_z_zz_y_x, g_y_0_x_0_z_zz_y_y, g_y_0_x_0_z_zz_y_z, g_yz_zz_xy_x, g_yz_zz_xy_y, g_yz_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_zz_y_x[i] = 4.0 * g_yz_zz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_y_y[i] = 4.0 * g_yz_zz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_y_z[i] = 4.0 * g_yz_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (645-648)

    #pragma omp simd aligned(g_y_0_x_0_z_zz_z_x, g_y_0_x_0_z_zz_z_y, g_y_0_x_0_z_zz_z_z, g_yz_zz_xz_x, g_yz_zz_xz_y, g_yz_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_zz_z_x[i] = 4.0 * g_yz_zz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_z_y[i] = 4.0 * g_yz_zz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_z_z[i] = 4.0 * g_yz_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (648-651)

    #pragma omp simd aligned(g_xy_xx_xy_x, g_xy_xx_xy_y, g_xy_xx_xy_z, g_y_0_y_0_x_xx_x_x, g_y_0_y_0_x_xx_x_y, g_y_0_y_0_x_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xx_x_x[i] = 4.0 * g_xy_xx_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_x_y[i] = 4.0 * g_xy_xx_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_x_z[i] = 4.0 * g_xy_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (651-654)

    #pragma omp simd aligned(g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xy_xx_yy_x, g_xy_xx_yy_y, g_xy_xx_yy_z, g_y_0_y_0_x_xx_y_x, g_y_0_y_0_x_xx_y_y, g_y_0_y_0_x_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xx_y_x[i] = -2.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xy_xx_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_y_y[i] = -2.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xy_xx_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_y_z[i] = -2.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xy_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (654-657)

    #pragma omp simd aligned(g_xy_xx_yz_x, g_xy_xx_yz_y, g_xy_xx_yz_z, g_y_0_y_0_x_xx_z_x, g_y_0_y_0_x_xx_z_y, g_y_0_y_0_x_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xx_z_x[i] = 4.0 * g_xy_xx_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_z_y[i] = 4.0 * g_xy_xx_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_z_z[i] = 4.0 * g_xy_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (657-660)

    #pragma omp simd aligned(g_xy_xy_xy_x, g_xy_xy_xy_y, g_xy_xy_xy_z, g_y_0_y_0_x_xy_x_x, g_y_0_y_0_x_xy_x_y, g_y_0_y_0_x_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xy_x_x[i] = 4.0 * g_xy_xy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_x_y[i] = 4.0 * g_xy_xy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_x_z[i] = 4.0 * g_xy_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (660-663)

    #pragma omp simd aligned(g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xy_xy_yy_x, g_xy_xy_yy_y, g_xy_xy_yy_z, g_y_0_y_0_x_xy_y_x, g_y_0_y_0_x_xy_y_y, g_y_0_y_0_x_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xy_y_x[i] = -2.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xy_xy_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_y_y[i] = -2.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xy_xy_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_y_z[i] = -2.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xy_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (663-666)

    #pragma omp simd aligned(g_xy_xy_yz_x, g_xy_xy_yz_y, g_xy_xy_yz_z, g_y_0_y_0_x_xy_z_x, g_y_0_y_0_x_xy_z_y, g_y_0_y_0_x_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xy_z_x[i] = 4.0 * g_xy_xy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_z_y[i] = 4.0 * g_xy_xy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_z_z[i] = 4.0 * g_xy_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (666-669)

    #pragma omp simd aligned(g_xy_xz_xy_x, g_xy_xz_xy_y, g_xy_xz_xy_z, g_y_0_y_0_x_xz_x_x, g_y_0_y_0_x_xz_x_y, g_y_0_y_0_x_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xz_x_x[i] = 4.0 * g_xy_xz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_x_y[i] = 4.0 * g_xy_xz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_x_z[i] = 4.0 * g_xy_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (669-672)

    #pragma omp simd aligned(g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xy_xz_yy_x, g_xy_xz_yy_y, g_xy_xz_yy_z, g_y_0_y_0_x_xz_y_x, g_y_0_y_0_x_xz_y_y, g_y_0_y_0_x_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xz_y_x[i] = -2.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xy_xz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_y_y[i] = -2.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xy_xz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_y_z[i] = -2.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xy_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (672-675)

    #pragma omp simd aligned(g_xy_xz_yz_x, g_xy_xz_yz_y, g_xy_xz_yz_z, g_y_0_y_0_x_xz_z_x, g_y_0_y_0_x_xz_z_y, g_y_0_y_0_x_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xz_z_x[i] = 4.0 * g_xy_xz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_z_y[i] = 4.0 * g_xy_xz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_z_z[i] = 4.0 * g_xy_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (675-678)

    #pragma omp simd aligned(g_xy_yy_xy_x, g_xy_yy_xy_y, g_xy_yy_xy_z, g_y_0_y_0_x_yy_x_x, g_y_0_y_0_x_yy_x_y, g_y_0_y_0_x_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_yy_x_x[i] = 4.0 * g_xy_yy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_x_y[i] = 4.0 * g_xy_yy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_x_z[i] = 4.0 * g_xy_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (678-681)

    #pragma omp simd aligned(g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xy_yy_yy_x, g_xy_yy_yy_y, g_xy_yy_yy_z, g_y_0_y_0_x_yy_y_x, g_y_0_y_0_x_yy_y_y, g_y_0_y_0_x_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_yy_y_x[i] = -2.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xy_yy_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_y_y[i] = -2.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xy_yy_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_y_z[i] = -2.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xy_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (681-684)

    #pragma omp simd aligned(g_xy_yy_yz_x, g_xy_yy_yz_y, g_xy_yy_yz_z, g_y_0_y_0_x_yy_z_x, g_y_0_y_0_x_yy_z_y, g_y_0_y_0_x_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_yy_z_x[i] = 4.0 * g_xy_yy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_z_y[i] = 4.0 * g_xy_yy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_z_z[i] = 4.0 * g_xy_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (684-687)

    #pragma omp simd aligned(g_xy_yz_xy_x, g_xy_yz_xy_y, g_xy_yz_xy_z, g_y_0_y_0_x_yz_x_x, g_y_0_y_0_x_yz_x_y, g_y_0_y_0_x_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_yz_x_x[i] = 4.0 * g_xy_yz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_x_y[i] = 4.0 * g_xy_yz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_x_z[i] = 4.0 * g_xy_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (687-690)

    #pragma omp simd aligned(g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xy_yz_yy_x, g_xy_yz_yy_y, g_xy_yz_yy_z, g_y_0_y_0_x_yz_y_x, g_y_0_y_0_x_yz_y_y, g_y_0_y_0_x_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_yz_y_x[i] = -2.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xy_yz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_y_y[i] = -2.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xy_yz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_y_z[i] = -2.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xy_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (690-693)

    #pragma omp simd aligned(g_xy_yz_yz_x, g_xy_yz_yz_y, g_xy_yz_yz_z, g_y_0_y_0_x_yz_z_x, g_y_0_y_0_x_yz_z_y, g_y_0_y_0_x_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_yz_z_x[i] = 4.0 * g_xy_yz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_z_y[i] = 4.0 * g_xy_yz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_z_z[i] = 4.0 * g_xy_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (693-696)

    #pragma omp simd aligned(g_xy_zz_xy_x, g_xy_zz_xy_y, g_xy_zz_xy_z, g_y_0_y_0_x_zz_x_x, g_y_0_y_0_x_zz_x_y, g_y_0_y_0_x_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_zz_x_x[i] = 4.0 * g_xy_zz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_x_y[i] = 4.0 * g_xy_zz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_x_z[i] = 4.0 * g_xy_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (696-699)

    #pragma omp simd aligned(g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xy_zz_yy_x, g_xy_zz_yy_y, g_xy_zz_yy_z, g_y_0_y_0_x_zz_y_x, g_y_0_y_0_x_zz_y_y, g_y_0_y_0_x_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_zz_y_x[i] = -2.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xy_zz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_y_y[i] = -2.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xy_zz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_y_z[i] = -2.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xy_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (699-702)

    #pragma omp simd aligned(g_xy_zz_yz_x, g_xy_zz_yz_y, g_xy_zz_yz_z, g_y_0_y_0_x_zz_z_x, g_y_0_y_0_x_zz_z_y, g_y_0_y_0_x_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_zz_z_x[i] = 4.0 * g_xy_zz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_z_y[i] = 4.0 * g_xy_zz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_z_z[i] = 4.0 * g_xy_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (702-705)

    #pragma omp simd aligned(g_0_xx_xy_x, g_0_xx_xy_y, g_0_xx_xy_z, g_y_0_y_0_y_xx_x_x, g_y_0_y_0_y_xx_x_y, g_y_0_y_0_y_xx_x_z, g_yy_xx_xy_x, g_yy_xx_xy_y, g_yy_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xx_x_x[i] = -2.0 * g_0_xx_xy_x[i] * c_exps[i] + 4.0 * g_yy_xx_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_x_y[i] = -2.0 * g_0_xx_xy_y[i] * c_exps[i] + 4.0 * g_yy_xx_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_x_z[i] = -2.0 * g_0_xx_xy_z[i] * c_exps[i] + 4.0 * g_yy_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (705-708)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_0_xx_yy_x, g_0_xx_yy_y, g_0_xx_yy_z, g_y_0_y_0_y_xx_y_x, g_y_0_y_0_y_xx_y_y, g_y_0_y_0_y_xx_y_z, g_yy_xx_0_x, g_yy_xx_0_y, g_yy_xx_0_z, g_yy_xx_yy_x, g_yy_xx_yy_y, g_yy_xx_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xx_y_x[i] = g_0_xx_0_x[i] - 2.0 * g_0_xx_yy_x[i] * c_exps[i] - 2.0 * g_yy_xx_0_x[i] * a_exp + 4.0 * g_yy_xx_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_y_y[i] = g_0_xx_0_y[i] - 2.0 * g_0_xx_yy_y[i] * c_exps[i] - 2.0 * g_yy_xx_0_y[i] * a_exp + 4.0 * g_yy_xx_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_y_z[i] = g_0_xx_0_z[i] - 2.0 * g_0_xx_yy_z[i] * c_exps[i] - 2.0 * g_yy_xx_0_z[i] * a_exp + 4.0 * g_yy_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (708-711)

    #pragma omp simd aligned(g_0_xx_yz_x, g_0_xx_yz_y, g_0_xx_yz_z, g_y_0_y_0_y_xx_z_x, g_y_0_y_0_y_xx_z_y, g_y_0_y_0_y_xx_z_z, g_yy_xx_yz_x, g_yy_xx_yz_y, g_yy_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xx_z_x[i] = -2.0 * g_0_xx_yz_x[i] * c_exps[i] + 4.0 * g_yy_xx_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_z_y[i] = -2.0 * g_0_xx_yz_y[i] * c_exps[i] + 4.0 * g_yy_xx_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_z_z[i] = -2.0 * g_0_xx_yz_z[i] * c_exps[i] + 4.0 * g_yy_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (711-714)

    #pragma omp simd aligned(g_0_xy_xy_x, g_0_xy_xy_y, g_0_xy_xy_z, g_y_0_y_0_y_xy_x_x, g_y_0_y_0_y_xy_x_y, g_y_0_y_0_y_xy_x_z, g_yy_xy_xy_x, g_yy_xy_xy_y, g_yy_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xy_x_x[i] = -2.0 * g_0_xy_xy_x[i] * c_exps[i] + 4.0 * g_yy_xy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_x_y[i] = -2.0 * g_0_xy_xy_y[i] * c_exps[i] + 4.0 * g_yy_xy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_x_z[i] = -2.0 * g_0_xy_xy_z[i] * c_exps[i] + 4.0 * g_yy_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (714-717)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_0_xy_yy_x, g_0_xy_yy_y, g_0_xy_yy_z, g_y_0_y_0_y_xy_y_x, g_y_0_y_0_y_xy_y_y, g_y_0_y_0_y_xy_y_z, g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z, g_yy_xy_yy_x, g_yy_xy_yy_y, g_yy_xy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xy_y_x[i] = g_0_xy_0_x[i] - 2.0 * g_0_xy_yy_x[i] * c_exps[i] - 2.0 * g_yy_xy_0_x[i] * a_exp + 4.0 * g_yy_xy_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_y_y[i] = g_0_xy_0_y[i] - 2.0 * g_0_xy_yy_y[i] * c_exps[i] - 2.0 * g_yy_xy_0_y[i] * a_exp + 4.0 * g_yy_xy_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_y_z[i] = g_0_xy_0_z[i] - 2.0 * g_0_xy_yy_z[i] * c_exps[i] - 2.0 * g_yy_xy_0_z[i] * a_exp + 4.0 * g_yy_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (717-720)

    #pragma omp simd aligned(g_0_xy_yz_x, g_0_xy_yz_y, g_0_xy_yz_z, g_y_0_y_0_y_xy_z_x, g_y_0_y_0_y_xy_z_y, g_y_0_y_0_y_xy_z_z, g_yy_xy_yz_x, g_yy_xy_yz_y, g_yy_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xy_z_x[i] = -2.0 * g_0_xy_yz_x[i] * c_exps[i] + 4.0 * g_yy_xy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_z_y[i] = -2.0 * g_0_xy_yz_y[i] * c_exps[i] + 4.0 * g_yy_xy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_z_z[i] = -2.0 * g_0_xy_yz_z[i] * c_exps[i] + 4.0 * g_yy_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (720-723)

    #pragma omp simd aligned(g_0_xz_xy_x, g_0_xz_xy_y, g_0_xz_xy_z, g_y_0_y_0_y_xz_x_x, g_y_0_y_0_y_xz_x_y, g_y_0_y_0_y_xz_x_z, g_yy_xz_xy_x, g_yy_xz_xy_y, g_yy_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xz_x_x[i] = -2.0 * g_0_xz_xy_x[i] * c_exps[i] + 4.0 * g_yy_xz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_x_y[i] = -2.0 * g_0_xz_xy_y[i] * c_exps[i] + 4.0 * g_yy_xz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_x_z[i] = -2.0 * g_0_xz_xy_z[i] * c_exps[i] + 4.0 * g_yy_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (723-726)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_0_xz_yy_x, g_0_xz_yy_y, g_0_xz_yy_z, g_y_0_y_0_y_xz_y_x, g_y_0_y_0_y_xz_y_y, g_y_0_y_0_y_xz_y_z, g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z, g_yy_xz_yy_x, g_yy_xz_yy_y, g_yy_xz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xz_y_x[i] = g_0_xz_0_x[i] - 2.0 * g_0_xz_yy_x[i] * c_exps[i] - 2.0 * g_yy_xz_0_x[i] * a_exp + 4.0 * g_yy_xz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_y_y[i] = g_0_xz_0_y[i] - 2.0 * g_0_xz_yy_y[i] * c_exps[i] - 2.0 * g_yy_xz_0_y[i] * a_exp + 4.0 * g_yy_xz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_y_z[i] = g_0_xz_0_z[i] - 2.0 * g_0_xz_yy_z[i] * c_exps[i] - 2.0 * g_yy_xz_0_z[i] * a_exp + 4.0 * g_yy_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (726-729)

    #pragma omp simd aligned(g_0_xz_yz_x, g_0_xz_yz_y, g_0_xz_yz_z, g_y_0_y_0_y_xz_z_x, g_y_0_y_0_y_xz_z_y, g_y_0_y_0_y_xz_z_z, g_yy_xz_yz_x, g_yy_xz_yz_y, g_yy_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xz_z_x[i] = -2.0 * g_0_xz_yz_x[i] * c_exps[i] + 4.0 * g_yy_xz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_z_y[i] = -2.0 * g_0_xz_yz_y[i] * c_exps[i] + 4.0 * g_yy_xz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_z_z[i] = -2.0 * g_0_xz_yz_z[i] * c_exps[i] + 4.0 * g_yy_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (729-732)

    #pragma omp simd aligned(g_0_yy_xy_x, g_0_yy_xy_y, g_0_yy_xy_z, g_y_0_y_0_y_yy_x_x, g_y_0_y_0_y_yy_x_y, g_y_0_y_0_y_yy_x_z, g_yy_yy_xy_x, g_yy_yy_xy_y, g_yy_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_yy_x_x[i] = -2.0 * g_0_yy_xy_x[i] * c_exps[i] + 4.0 * g_yy_yy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_x_y[i] = -2.0 * g_0_yy_xy_y[i] * c_exps[i] + 4.0 * g_yy_yy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_x_z[i] = -2.0 * g_0_yy_xy_z[i] * c_exps[i] + 4.0 * g_yy_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (732-735)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_0_yy_yy_x, g_0_yy_yy_y, g_0_yy_yy_z, g_y_0_y_0_y_yy_y_x, g_y_0_y_0_y_yy_y_y, g_y_0_y_0_y_yy_y_z, g_yy_yy_0_x, g_yy_yy_0_y, g_yy_yy_0_z, g_yy_yy_yy_x, g_yy_yy_yy_y, g_yy_yy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_yy_y_x[i] = g_0_yy_0_x[i] - 2.0 * g_0_yy_yy_x[i] * c_exps[i] - 2.0 * g_yy_yy_0_x[i] * a_exp + 4.0 * g_yy_yy_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_y_y[i] = g_0_yy_0_y[i] - 2.0 * g_0_yy_yy_y[i] * c_exps[i] - 2.0 * g_yy_yy_0_y[i] * a_exp + 4.0 * g_yy_yy_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_y_z[i] = g_0_yy_0_z[i] - 2.0 * g_0_yy_yy_z[i] * c_exps[i] - 2.0 * g_yy_yy_0_z[i] * a_exp + 4.0 * g_yy_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (735-738)

    #pragma omp simd aligned(g_0_yy_yz_x, g_0_yy_yz_y, g_0_yy_yz_z, g_y_0_y_0_y_yy_z_x, g_y_0_y_0_y_yy_z_y, g_y_0_y_0_y_yy_z_z, g_yy_yy_yz_x, g_yy_yy_yz_y, g_yy_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_yy_z_x[i] = -2.0 * g_0_yy_yz_x[i] * c_exps[i] + 4.0 * g_yy_yy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_z_y[i] = -2.0 * g_0_yy_yz_y[i] * c_exps[i] + 4.0 * g_yy_yy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_z_z[i] = -2.0 * g_0_yy_yz_z[i] * c_exps[i] + 4.0 * g_yy_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (738-741)

    #pragma omp simd aligned(g_0_yz_xy_x, g_0_yz_xy_y, g_0_yz_xy_z, g_y_0_y_0_y_yz_x_x, g_y_0_y_0_y_yz_x_y, g_y_0_y_0_y_yz_x_z, g_yy_yz_xy_x, g_yy_yz_xy_y, g_yy_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_yz_x_x[i] = -2.0 * g_0_yz_xy_x[i] * c_exps[i] + 4.0 * g_yy_yz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_x_y[i] = -2.0 * g_0_yz_xy_y[i] * c_exps[i] + 4.0 * g_yy_yz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_x_z[i] = -2.0 * g_0_yz_xy_z[i] * c_exps[i] + 4.0 * g_yy_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (741-744)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_0_yz_yy_x, g_0_yz_yy_y, g_0_yz_yy_z, g_y_0_y_0_y_yz_y_x, g_y_0_y_0_y_yz_y_y, g_y_0_y_0_y_yz_y_z, g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z, g_yy_yz_yy_x, g_yy_yz_yy_y, g_yy_yz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_yz_y_x[i] = g_0_yz_0_x[i] - 2.0 * g_0_yz_yy_x[i] * c_exps[i] - 2.0 * g_yy_yz_0_x[i] * a_exp + 4.0 * g_yy_yz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_y_y[i] = g_0_yz_0_y[i] - 2.0 * g_0_yz_yy_y[i] * c_exps[i] - 2.0 * g_yy_yz_0_y[i] * a_exp + 4.0 * g_yy_yz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_y_z[i] = g_0_yz_0_z[i] - 2.0 * g_0_yz_yy_z[i] * c_exps[i] - 2.0 * g_yy_yz_0_z[i] * a_exp + 4.0 * g_yy_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (744-747)

    #pragma omp simd aligned(g_0_yz_yz_x, g_0_yz_yz_y, g_0_yz_yz_z, g_y_0_y_0_y_yz_z_x, g_y_0_y_0_y_yz_z_y, g_y_0_y_0_y_yz_z_z, g_yy_yz_yz_x, g_yy_yz_yz_y, g_yy_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_yz_z_x[i] = -2.0 * g_0_yz_yz_x[i] * c_exps[i] + 4.0 * g_yy_yz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_z_y[i] = -2.0 * g_0_yz_yz_y[i] * c_exps[i] + 4.0 * g_yy_yz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_z_z[i] = -2.0 * g_0_yz_yz_z[i] * c_exps[i] + 4.0 * g_yy_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (747-750)

    #pragma omp simd aligned(g_0_zz_xy_x, g_0_zz_xy_y, g_0_zz_xy_z, g_y_0_y_0_y_zz_x_x, g_y_0_y_0_y_zz_x_y, g_y_0_y_0_y_zz_x_z, g_yy_zz_xy_x, g_yy_zz_xy_y, g_yy_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_zz_x_x[i] = -2.0 * g_0_zz_xy_x[i] * c_exps[i] + 4.0 * g_yy_zz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_x_y[i] = -2.0 * g_0_zz_xy_y[i] * c_exps[i] + 4.0 * g_yy_zz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_x_z[i] = -2.0 * g_0_zz_xy_z[i] * c_exps[i] + 4.0 * g_yy_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (750-753)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_0_zz_yy_x, g_0_zz_yy_y, g_0_zz_yy_z, g_y_0_y_0_y_zz_y_x, g_y_0_y_0_y_zz_y_y, g_y_0_y_0_y_zz_y_z, g_yy_zz_0_x, g_yy_zz_0_y, g_yy_zz_0_z, g_yy_zz_yy_x, g_yy_zz_yy_y, g_yy_zz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_zz_y_x[i] = g_0_zz_0_x[i] - 2.0 * g_0_zz_yy_x[i] * c_exps[i] - 2.0 * g_yy_zz_0_x[i] * a_exp + 4.0 * g_yy_zz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_y_y[i] = g_0_zz_0_y[i] - 2.0 * g_0_zz_yy_y[i] * c_exps[i] - 2.0 * g_yy_zz_0_y[i] * a_exp + 4.0 * g_yy_zz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_y_z[i] = g_0_zz_0_z[i] - 2.0 * g_0_zz_yy_z[i] * c_exps[i] - 2.0 * g_yy_zz_0_z[i] * a_exp + 4.0 * g_yy_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (753-756)

    #pragma omp simd aligned(g_0_zz_yz_x, g_0_zz_yz_y, g_0_zz_yz_z, g_y_0_y_0_y_zz_z_x, g_y_0_y_0_y_zz_z_y, g_y_0_y_0_y_zz_z_z, g_yy_zz_yz_x, g_yy_zz_yz_y, g_yy_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_zz_z_x[i] = -2.0 * g_0_zz_yz_x[i] * c_exps[i] + 4.0 * g_yy_zz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_z_y[i] = -2.0 * g_0_zz_yz_y[i] * c_exps[i] + 4.0 * g_yy_zz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_z_z[i] = -2.0 * g_0_zz_yz_z[i] * c_exps[i] + 4.0 * g_yy_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (756-759)

    #pragma omp simd aligned(g_y_0_y_0_z_xx_x_x, g_y_0_y_0_z_xx_x_y, g_y_0_y_0_z_xx_x_z, g_yz_xx_xy_x, g_yz_xx_xy_y, g_yz_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xx_x_x[i] = 4.0 * g_yz_xx_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_x_y[i] = 4.0 * g_yz_xx_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_x_z[i] = 4.0 * g_yz_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (759-762)

    #pragma omp simd aligned(g_y_0_y_0_z_xx_y_x, g_y_0_y_0_z_xx_y_y, g_y_0_y_0_z_xx_y_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z, g_yz_xx_yy_x, g_yz_xx_yy_y, g_yz_xx_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xx_y_x[i] = -2.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_yz_xx_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_y_y[i] = -2.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_yz_xx_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_y_z[i] = -2.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_yz_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (762-765)

    #pragma omp simd aligned(g_y_0_y_0_z_xx_z_x, g_y_0_y_0_z_xx_z_y, g_y_0_y_0_z_xx_z_z, g_yz_xx_yz_x, g_yz_xx_yz_y, g_yz_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xx_z_x[i] = 4.0 * g_yz_xx_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_z_y[i] = 4.0 * g_yz_xx_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_z_z[i] = 4.0 * g_yz_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (765-768)

    #pragma omp simd aligned(g_y_0_y_0_z_xy_x_x, g_y_0_y_0_z_xy_x_y, g_y_0_y_0_z_xy_x_z, g_yz_xy_xy_x, g_yz_xy_xy_y, g_yz_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xy_x_x[i] = 4.0 * g_yz_xy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_x_y[i] = 4.0 * g_yz_xy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_x_z[i] = 4.0 * g_yz_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (768-771)

    #pragma omp simd aligned(g_y_0_y_0_z_xy_y_x, g_y_0_y_0_z_xy_y_y, g_y_0_y_0_z_xy_y_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_yz_xy_yy_x, g_yz_xy_yy_y, g_yz_xy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xy_y_x[i] = -2.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_yz_xy_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_y_y[i] = -2.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_yz_xy_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_y_z[i] = -2.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_yz_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (771-774)

    #pragma omp simd aligned(g_y_0_y_0_z_xy_z_x, g_y_0_y_0_z_xy_z_y, g_y_0_y_0_z_xy_z_z, g_yz_xy_yz_x, g_yz_xy_yz_y, g_yz_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xy_z_x[i] = 4.0 * g_yz_xy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_z_y[i] = 4.0 * g_yz_xy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_z_z[i] = 4.0 * g_yz_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (774-777)

    #pragma omp simd aligned(g_y_0_y_0_z_xz_x_x, g_y_0_y_0_z_xz_x_y, g_y_0_y_0_z_xz_x_z, g_yz_xz_xy_x, g_yz_xz_xy_y, g_yz_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xz_x_x[i] = 4.0 * g_yz_xz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_x_y[i] = 4.0 * g_yz_xz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_x_z[i] = 4.0 * g_yz_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (777-780)

    #pragma omp simd aligned(g_y_0_y_0_z_xz_y_x, g_y_0_y_0_z_xz_y_y, g_y_0_y_0_z_xz_y_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_yz_xz_yy_x, g_yz_xz_yy_y, g_yz_xz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xz_y_x[i] = -2.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_yz_xz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_y_y[i] = -2.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_yz_xz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_y_z[i] = -2.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_yz_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (780-783)

    #pragma omp simd aligned(g_y_0_y_0_z_xz_z_x, g_y_0_y_0_z_xz_z_y, g_y_0_y_0_z_xz_z_z, g_yz_xz_yz_x, g_yz_xz_yz_y, g_yz_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xz_z_x[i] = 4.0 * g_yz_xz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_z_y[i] = 4.0 * g_yz_xz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_z_z[i] = 4.0 * g_yz_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (783-786)

    #pragma omp simd aligned(g_y_0_y_0_z_yy_x_x, g_y_0_y_0_z_yy_x_y, g_y_0_y_0_z_yy_x_z, g_yz_yy_xy_x, g_yz_yy_xy_y, g_yz_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_yy_x_x[i] = 4.0 * g_yz_yy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_x_y[i] = 4.0 * g_yz_yy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_x_z[i] = 4.0 * g_yz_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (786-789)

    #pragma omp simd aligned(g_y_0_y_0_z_yy_y_x, g_y_0_y_0_z_yy_y_y, g_y_0_y_0_z_yy_y_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z, g_yz_yy_yy_x, g_yz_yy_yy_y, g_yz_yy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_yy_y_x[i] = -2.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_yz_yy_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_y_y[i] = -2.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_yz_yy_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_y_z[i] = -2.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_yz_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (789-792)

    #pragma omp simd aligned(g_y_0_y_0_z_yy_z_x, g_y_0_y_0_z_yy_z_y, g_y_0_y_0_z_yy_z_z, g_yz_yy_yz_x, g_yz_yy_yz_y, g_yz_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_yy_z_x[i] = 4.0 * g_yz_yy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_z_y[i] = 4.0 * g_yz_yy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_z_z[i] = 4.0 * g_yz_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (792-795)

    #pragma omp simd aligned(g_y_0_y_0_z_yz_x_x, g_y_0_y_0_z_yz_x_y, g_y_0_y_0_z_yz_x_z, g_yz_yz_xy_x, g_yz_yz_xy_y, g_yz_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_yz_x_x[i] = 4.0 * g_yz_yz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_x_y[i] = 4.0 * g_yz_yz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_x_z[i] = 4.0 * g_yz_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (795-798)

    #pragma omp simd aligned(g_y_0_y_0_z_yz_y_x, g_y_0_y_0_z_yz_y_y, g_y_0_y_0_z_yz_y_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_yz_yz_yy_x, g_yz_yz_yy_y, g_yz_yz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_yz_y_x[i] = -2.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_yz_yz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_y_y[i] = -2.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_yz_yz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_y_z[i] = -2.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_yz_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (798-801)

    #pragma omp simd aligned(g_y_0_y_0_z_yz_z_x, g_y_0_y_0_z_yz_z_y, g_y_0_y_0_z_yz_z_z, g_yz_yz_yz_x, g_yz_yz_yz_y, g_yz_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_yz_z_x[i] = 4.0 * g_yz_yz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_z_y[i] = 4.0 * g_yz_yz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_z_z[i] = 4.0 * g_yz_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (801-804)

    #pragma omp simd aligned(g_y_0_y_0_z_zz_x_x, g_y_0_y_0_z_zz_x_y, g_y_0_y_0_z_zz_x_z, g_yz_zz_xy_x, g_yz_zz_xy_y, g_yz_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_zz_x_x[i] = 4.0 * g_yz_zz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_x_y[i] = 4.0 * g_yz_zz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_x_z[i] = 4.0 * g_yz_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (804-807)

    #pragma omp simd aligned(g_y_0_y_0_z_zz_y_x, g_y_0_y_0_z_zz_y_y, g_y_0_y_0_z_zz_y_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z, g_yz_zz_yy_x, g_yz_zz_yy_y, g_yz_zz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_zz_y_x[i] = -2.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_yz_zz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_y_y[i] = -2.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_yz_zz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_y_z[i] = -2.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_yz_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (807-810)

    #pragma omp simd aligned(g_y_0_y_0_z_zz_z_x, g_y_0_y_0_z_zz_z_y, g_y_0_y_0_z_zz_z_z, g_yz_zz_yz_x, g_yz_zz_yz_y, g_yz_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_zz_z_x[i] = 4.0 * g_yz_zz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_z_y[i] = 4.0 * g_yz_zz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_z_z[i] = 4.0 * g_yz_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (810-813)

    #pragma omp simd aligned(g_xy_xx_xz_x, g_xy_xx_xz_y, g_xy_xx_xz_z, g_y_0_z_0_x_xx_x_x, g_y_0_z_0_x_xx_x_y, g_y_0_z_0_x_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xx_x_x[i] = 4.0 * g_xy_xx_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_x_y[i] = 4.0 * g_xy_xx_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_x_z[i] = 4.0 * g_xy_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (813-816)

    #pragma omp simd aligned(g_xy_xx_yz_x, g_xy_xx_yz_y, g_xy_xx_yz_z, g_y_0_z_0_x_xx_y_x, g_y_0_z_0_x_xx_y_y, g_y_0_z_0_x_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xx_y_x[i] = 4.0 * g_xy_xx_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_y_y[i] = 4.0 * g_xy_xx_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_y_z[i] = 4.0 * g_xy_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (816-819)

    #pragma omp simd aligned(g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xy_xx_zz_x, g_xy_xx_zz_y, g_xy_xx_zz_z, g_y_0_z_0_x_xx_z_x, g_y_0_z_0_x_xx_z_y, g_y_0_z_0_x_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xx_z_x[i] = -2.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xy_xx_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_z_y[i] = -2.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xy_xx_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_z_z[i] = -2.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xy_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (819-822)

    #pragma omp simd aligned(g_xy_xy_xz_x, g_xy_xy_xz_y, g_xy_xy_xz_z, g_y_0_z_0_x_xy_x_x, g_y_0_z_0_x_xy_x_y, g_y_0_z_0_x_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xy_x_x[i] = 4.0 * g_xy_xy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_x_y[i] = 4.0 * g_xy_xy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_x_z[i] = 4.0 * g_xy_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (822-825)

    #pragma omp simd aligned(g_xy_xy_yz_x, g_xy_xy_yz_y, g_xy_xy_yz_z, g_y_0_z_0_x_xy_y_x, g_y_0_z_0_x_xy_y_y, g_y_0_z_0_x_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xy_y_x[i] = 4.0 * g_xy_xy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_y_y[i] = 4.0 * g_xy_xy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_y_z[i] = 4.0 * g_xy_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (825-828)

    #pragma omp simd aligned(g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xy_xy_zz_x, g_xy_xy_zz_y, g_xy_xy_zz_z, g_y_0_z_0_x_xy_z_x, g_y_0_z_0_x_xy_z_y, g_y_0_z_0_x_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xy_z_x[i] = -2.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xy_xy_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_z_y[i] = -2.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xy_xy_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_z_z[i] = -2.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xy_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (828-831)

    #pragma omp simd aligned(g_xy_xz_xz_x, g_xy_xz_xz_y, g_xy_xz_xz_z, g_y_0_z_0_x_xz_x_x, g_y_0_z_0_x_xz_x_y, g_y_0_z_0_x_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xz_x_x[i] = 4.0 * g_xy_xz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_x_y[i] = 4.0 * g_xy_xz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_x_z[i] = 4.0 * g_xy_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (831-834)

    #pragma omp simd aligned(g_xy_xz_yz_x, g_xy_xz_yz_y, g_xy_xz_yz_z, g_y_0_z_0_x_xz_y_x, g_y_0_z_0_x_xz_y_y, g_y_0_z_0_x_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xz_y_x[i] = 4.0 * g_xy_xz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_y_y[i] = 4.0 * g_xy_xz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_y_z[i] = 4.0 * g_xy_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (834-837)

    #pragma omp simd aligned(g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xy_xz_zz_x, g_xy_xz_zz_y, g_xy_xz_zz_z, g_y_0_z_0_x_xz_z_x, g_y_0_z_0_x_xz_z_y, g_y_0_z_0_x_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xz_z_x[i] = -2.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xy_xz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_z_y[i] = -2.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xy_xz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_z_z[i] = -2.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xy_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (837-840)

    #pragma omp simd aligned(g_xy_yy_xz_x, g_xy_yy_xz_y, g_xy_yy_xz_z, g_y_0_z_0_x_yy_x_x, g_y_0_z_0_x_yy_x_y, g_y_0_z_0_x_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_yy_x_x[i] = 4.0 * g_xy_yy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_x_y[i] = 4.0 * g_xy_yy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_x_z[i] = 4.0 * g_xy_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (840-843)

    #pragma omp simd aligned(g_xy_yy_yz_x, g_xy_yy_yz_y, g_xy_yy_yz_z, g_y_0_z_0_x_yy_y_x, g_y_0_z_0_x_yy_y_y, g_y_0_z_0_x_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_yy_y_x[i] = 4.0 * g_xy_yy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_y_y[i] = 4.0 * g_xy_yy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_y_z[i] = 4.0 * g_xy_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (843-846)

    #pragma omp simd aligned(g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xy_yy_zz_x, g_xy_yy_zz_y, g_xy_yy_zz_z, g_y_0_z_0_x_yy_z_x, g_y_0_z_0_x_yy_z_y, g_y_0_z_0_x_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_yy_z_x[i] = -2.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xy_yy_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_z_y[i] = -2.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xy_yy_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_z_z[i] = -2.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xy_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (846-849)

    #pragma omp simd aligned(g_xy_yz_xz_x, g_xy_yz_xz_y, g_xy_yz_xz_z, g_y_0_z_0_x_yz_x_x, g_y_0_z_0_x_yz_x_y, g_y_0_z_0_x_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_yz_x_x[i] = 4.0 * g_xy_yz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_x_y[i] = 4.0 * g_xy_yz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_x_z[i] = 4.0 * g_xy_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (849-852)

    #pragma omp simd aligned(g_xy_yz_yz_x, g_xy_yz_yz_y, g_xy_yz_yz_z, g_y_0_z_0_x_yz_y_x, g_y_0_z_0_x_yz_y_y, g_y_0_z_0_x_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_yz_y_x[i] = 4.0 * g_xy_yz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_y_y[i] = 4.0 * g_xy_yz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_y_z[i] = 4.0 * g_xy_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (852-855)

    #pragma omp simd aligned(g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xy_yz_zz_x, g_xy_yz_zz_y, g_xy_yz_zz_z, g_y_0_z_0_x_yz_z_x, g_y_0_z_0_x_yz_z_y, g_y_0_z_0_x_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_yz_z_x[i] = -2.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xy_yz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_z_y[i] = -2.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xy_yz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_z_z[i] = -2.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xy_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (855-858)

    #pragma omp simd aligned(g_xy_zz_xz_x, g_xy_zz_xz_y, g_xy_zz_xz_z, g_y_0_z_0_x_zz_x_x, g_y_0_z_0_x_zz_x_y, g_y_0_z_0_x_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_zz_x_x[i] = 4.0 * g_xy_zz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_x_y[i] = 4.0 * g_xy_zz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_x_z[i] = 4.0 * g_xy_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (858-861)

    #pragma omp simd aligned(g_xy_zz_yz_x, g_xy_zz_yz_y, g_xy_zz_yz_z, g_y_0_z_0_x_zz_y_x, g_y_0_z_0_x_zz_y_y, g_y_0_z_0_x_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_zz_y_x[i] = 4.0 * g_xy_zz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_y_y[i] = 4.0 * g_xy_zz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_y_z[i] = 4.0 * g_xy_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (861-864)

    #pragma omp simd aligned(g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xy_zz_zz_x, g_xy_zz_zz_y, g_xy_zz_zz_z, g_y_0_z_0_x_zz_z_x, g_y_0_z_0_x_zz_z_y, g_y_0_z_0_x_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_zz_z_x[i] = -2.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xy_zz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_z_y[i] = -2.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xy_zz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_z_z[i] = -2.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xy_zz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (864-867)

    #pragma omp simd aligned(g_0_xx_xz_x, g_0_xx_xz_y, g_0_xx_xz_z, g_y_0_z_0_y_xx_x_x, g_y_0_z_0_y_xx_x_y, g_y_0_z_0_y_xx_x_z, g_yy_xx_xz_x, g_yy_xx_xz_y, g_yy_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xx_x_x[i] = -2.0 * g_0_xx_xz_x[i] * c_exps[i] + 4.0 * g_yy_xx_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_x_y[i] = -2.0 * g_0_xx_xz_y[i] * c_exps[i] + 4.0 * g_yy_xx_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_x_z[i] = -2.0 * g_0_xx_xz_z[i] * c_exps[i] + 4.0 * g_yy_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (867-870)

    #pragma omp simd aligned(g_0_xx_yz_x, g_0_xx_yz_y, g_0_xx_yz_z, g_y_0_z_0_y_xx_y_x, g_y_0_z_0_y_xx_y_y, g_y_0_z_0_y_xx_y_z, g_yy_xx_yz_x, g_yy_xx_yz_y, g_yy_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xx_y_x[i] = -2.0 * g_0_xx_yz_x[i] * c_exps[i] + 4.0 * g_yy_xx_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_y_y[i] = -2.0 * g_0_xx_yz_y[i] * c_exps[i] + 4.0 * g_yy_xx_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_y_z[i] = -2.0 * g_0_xx_yz_z[i] * c_exps[i] + 4.0 * g_yy_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (870-873)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_0_xx_zz_x, g_0_xx_zz_y, g_0_xx_zz_z, g_y_0_z_0_y_xx_z_x, g_y_0_z_0_y_xx_z_y, g_y_0_z_0_y_xx_z_z, g_yy_xx_0_x, g_yy_xx_0_y, g_yy_xx_0_z, g_yy_xx_zz_x, g_yy_xx_zz_y, g_yy_xx_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xx_z_x[i] = g_0_xx_0_x[i] - 2.0 * g_0_xx_zz_x[i] * c_exps[i] - 2.0 * g_yy_xx_0_x[i] * a_exp + 4.0 * g_yy_xx_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_z_y[i] = g_0_xx_0_y[i] - 2.0 * g_0_xx_zz_y[i] * c_exps[i] - 2.0 * g_yy_xx_0_y[i] * a_exp + 4.0 * g_yy_xx_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_z_z[i] = g_0_xx_0_z[i] - 2.0 * g_0_xx_zz_z[i] * c_exps[i] - 2.0 * g_yy_xx_0_z[i] * a_exp + 4.0 * g_yy_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (873-876)

    #pragma omp simd aligned(g_0_xy_xz_x, g_0_xy_xz_y, g_0_xy_xz_z, g_y_0_z_0_y_xy_x_x, g_y_0_z_0_y_xy_x_y, g_y_0_z_0_y_xy_x_z, g_yy_xy_xz_x, g_yy_xy_xz_y, g_yy_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xy_x_x[i] = -2.0 * g_0_xy_xz_x[i] * c_exps[i] + 4.0 * g_yy_xy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_x_y[i] = -2.0 * g_0_xy_xz_y[i] * c_exps[i] + 4.0 * g_yy_xy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_x_z[i] = -2.0 * g_0_xy_xz_z[i] * c_exps[i] + 4.0 * g_yy_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (876-879)

    #pragma omp simd aligned(g_0_xy_yz_x, g_0_xy_yz_y, g_0_xy_yz_z, g_y_0_z_0_y_xy_y_x, g_y_0_z_0_y_xy_y_y, g_y_0_z_0_y_xy_y_z, g_yy_xy_yz_x, g_yy_xy_yz_y, g_yy_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xy_y_x[i] = -2.0 * g_0_xy_yz_x[i] * c_exps[i] + 4.0 * g_yy_xy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_y_y[i] = -2.0 * g_0_xy_yz_y[i] * c_exps[i] + 4.0 * g_yy_xy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_y_z[i] = -2.0 * g_0_xy_yz_z[i] * c_exps[i] + 4.0 * g_yy_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (879-882)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_0_xy_zz_x, g_0_xy_zz_y, g_0_xy_zz_z, g_y_0_z_0_y_xy_z_x, g_y_0_z_0_y_xy_z_y, g_y_0_z_0_y_xy_z_z, g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z, g_yy_xy_zz_x, g_yy_xy_zz_y, g_yy_xy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xy_z_x[i] = g_0_xy_0_x[i] - 2.0 * g_0_xy_zz_x[i] * c_exps[i] - 2.0 * g_yy_xy_0_x[i] * a_exp + 4.0 * g_yy_xy_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_z_y[i] = g_0_xy_0_y[i] - 2.0 * g_0_xy_zz_y[i] * c_exps[i] - 2.0 * g_yy_xy_0_y[i] * a_exp + 4.0 * g_yy_xy_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_z_z[i] = g_0_xy_0_z[i] - 2.0 * g_0_xy_zz_z[i] * c_exps[i] - 2.0 * g_yy_xy_0_z[i] * a_exp + 4.0 * g_yy_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (882-885)

    #pragma omp simd aligned(g_0_xz_xz_x, g_0_xz_xz_y, g_0_xz_xz_z, g_y_0_z_0_y_xz_x_x, g_y_0_z_0_y_xz_x_y, g_y_0_z_0_y_xz_x_z, g_yy_xz_xz_x, g_yy_xz_xz_y, g_yy_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xz_x_x[i] = -2.0 * g_0_xz_xz_x[i] * c_exps[i] + 4.0 * g_yy_xz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_x_y[i] = -2.0 * g_0_xz_xz_y[i] * c_exps[i] + 4.0 * g_yy_xz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_x_z[i] = -2.0 * g_0_xz_xz_z[i] * c_exps[i] + 4.0 * g_yy_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (885-888)

    #pragma omp simd aligned(g_0_xz_yz_x, g_0_xz_yz_y, g_0_xz_yz_z, g_y_0_z_0_y_xz_y_x, g_y_0_z_0_y_xz_y_y, g_y_0_z_0_y_xz_y_z, g_yy_xz_yz_x, g_yy_xz_yz_y, g_yy_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xz_y_x[i] = -2.0 * g_0_xz_yz_x[i] * c_exps[i] + 4.0 * g_yy_xz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_y_y[i] = -2.0 * g_0_xz_yz_y[i] * c_exps[i] + 4.0 * g_yy_xz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_y_z[i] = -2.0 * g_0_xz_yz_z[i] * c_exps[i] + 4.0 * g_yy_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (888-891)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_0_xz_zz_x, g_0_xz_zz_y, g_0_xz_zz_z, g_y_0_z_0_y_xz_z_x, g_y_0_z_0_y_xz_z_y, g_y_0_z_0_y_xz_z_z, g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z, g_yy_xz_zz_x, g_yy_xz_zz_y, g_yy_xz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xz_z_x[i] = g_0_xz_0_x[i] - 2.0 * g_0_xz_zz_x[i] * c_exps[i] - 2.0 * g_yy_xz_0_x[i] * a_exp + 4.0 * g_yy_xz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_z_y[i] = g_0_xz_0_y[i] - 2.0 * g_0_xz_zz_y[i] * c_exps[i] - 2.0 * g_yy_xz_0_y[i] * a_exp + 4.0 * g_yy_xz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_z_z[i] = g_0_xz_0_z[i] - 2.0 * g_0_xz_zz_z[i] * c_exps[i] - 2.0 * g_yy_xz_0_z[i] * a_exp + 4.0 * g_yy_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (891-894)

    #pragma omp simd aligned(g_0_yy_xz_x, g_0_yy_xz_y, g_0_yy_xz_z, g_y_0_z_0_y_yy_x_x, g_y_0_z_0_y_yy_x_y, g_y_0_z_0_y_yy_x_z, g_yy_yy_xz_x, g_yy_yy_xz_y, g_yy_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_yy_x_x[i] = -2.0 * g_0_yy_xz_x[i] * c_exps[i] + 4.0 * g_yy_yy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_x_y[i] = -2.0 * g_0_yy_xz_y[i] * c_exps[i] + 4.0 * g_yy_yy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_x_z[i] = -2.0 * g_0_yy_xz_z[i] * c_exps[i] + 4.0 * g_yy_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (894-897)

    #pragma omp simd aligned(g_0_yy_yz_x, g_0_yy_yz_y, g_0_yy_yz_z, g_y_0_z_0_y_yy_y_x, g_y_0_z_0_y_yy_y_y, g_y_0_z_0_y_yy_y_z, g_yy_yy_yz_x, g_yy_yy_yz_y, g_yy_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_yy_y_x[i] = -2.0 * g_0_yy_yz_x[i] * c_exps[i] + 4.0 * g_yy_yy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_y_y[i] = -2.0 * g_0_yy_yz_y[i] * c_exps[i] + 4.0 * g_yy_yy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_y_z[i] = -2.0 * g_0_yy_yz_z[i] * c_exps[i] + 4.0 * g_yy_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (897-900)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_0_yy_zz_x, g_0_yy_zz_y, g_0_yy_zz_z, g_y_0_z_0_y_yy_z_x, g_y_0_z_0_y_yy_z_y, g_y_0_z_0_y_yy_z_z, g_yy_yy_0_x, g_yy_yy_0_y, g_yy_yy_0_z, g_yy_yy_zz_x, g_yy_yy_zz_y, g_yy_yy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_yy_z_x[i] = g_0_yy_0_x[i] - 2.0 * g_0_yy_zz_x[i] * c_exps[i] - 2.0 * g_yy_yy_0_x[i] * a_exp + 4.0 * g_yy_yy_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_z_y[i] = g_0_yy_0_y[i] - 2.0 * g_0_yy_zz_y[i] * c_exps[i] - 2.0 * g_yy_yy_0_y[i] * a_exp + 4.0 * g_yy_yy_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_z_z[i] = g_0_yy_0_z[i] - 2.0 * g_0_yy_zz_z[i] * c_exps[i] - 2.0 * g_yy_yy_0_z[i] * a_exp + 4.0 * g_yy_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (900-903)

    #pragma omp simd aligned(g_0_yz_xz_x, g_0_yz_xz_y, g_0_yz_xz_z, g_y_0_z_0_y_yz_x_x, g_y_0_z_0_y_yz_x_y, g_y_0_z_0_y_yz_x_z, g_yy_yz_xz_x, g_yy_yz_xz_y, g_yy_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_yz_x_x[i] = -2.0 * g_0_yz_xz_x[i] * c_exps[i] + 4.0 * g_yy_yz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_x_y[i] = -2.0 * g_0_yz_xz_y[i] * c_exps[i] + 4.0 * g_yy_yz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_x_z[i] = -2.0 * g_0_yz_xz_z[i] * c_exps[i] + 4.0 * g_yy_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (903-906)

    #pragma omp simd aligned(g_0_yz_yz_x, g_0_yz_yz_y, g_0_yz_yz_z, g_y_0_z_0_y_yz_y_x, g_y_0_z_0_y_yz_y_y, g_y_0_z_0_y_yz_y_z, g_yy_yz_yz_x, g_yy_yz_yz_y, g_yy_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_yz_y_x[i] = -2.0 * g_0_yz_yz_x[i] * c_exps[i] + 4.0 * g_yy_yz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_y_y[i] = -2.0 * g_0_yz_yz_y[i] * c_exps[i] + 4.0 * g_yy_yz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_y_z[i] = -2.0 * g_0_yz_yz_z[i] * c_exps[i] + 4.0 * g_yy_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (906-909)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_0_yz_zz_x, g_0_yz_zz_y, g_0_yz_zz_z, g_y_0_z_0_y_yz_z_x, g_y_0_z_0_y_yz_z_y, g_y_0_z_0_y_yz_z_z, g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z, g_yy_yz_zz_x, g_yy_yz_zz_y, g_yy_yz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_yz_z_x[i] = g_0_yz_0_x[i] - 2.0 * g_0_yz_zz_x[i] * c_exps[i] - 2.0 * g_yy_yz_0_x[i] * a_exp + 4.0 * g_yy_yz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_z_y[i] = g_0_yz_0_y[i] - 2.0 * g_0_yz_zz_y[i] * c_exps[i] - 2.0 * g_yy_yz_0_y[i] * a_exp + 4.0 * g_yy_yz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_z_z[i] = g_0_yz_0_z[i] - 2.0 * g_0_yz_zz_z[i] * c_exps[i] - 2.0 * g_yy_yz_0_z[i] * a_exp + 4.0 * g_yy_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (909-912)

    #pragma omp simd aligned(g_0_zz_xz_x, g_0_zz_xz_y, g_0_zz_xz_z, g_y_0_z_0_y_zz_x_x, g_y_0_z_0_y_zz_x_y, g_y_0_z_0_y_zz_x_z, g_yy_zz_xz_x, g_yy_zz_xz_y, g_yy_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_zz_x_x[i] = -2.0 * g_0_zz_xz_x[i] * c_exps[i] + 4.0 * g_yy_zz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_x_y[i] = -2.0 * g_0_zz_xz_y[i] * c_exps[i] + 4.0 * g_yy_zz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_x_z[i] = -2.0 * g_0_zz_xz_z[i] * c_exps[i] + 4.0 * g_yy_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (912-915)

    #pragma omp simd aligned(g_0_zz_yz_x, g_0_zz_yz_y, g_0_zz_yz_z, g_y_0_z_0_y_zz_y_x, g_y_0_z_0_y_zz_y_y, g_y_0_z_0_y_zz_y_z, g_yy_zz_yz_x, g_yy_zz_yz_y, g_yy_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_zz_y_x[i] = -2.0 * g_0_zz_yz_x[i] * c_exps[i] + 4.0 * g_yy_zz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_y_y[i] = -2.0 * g_0_zz_yz_y[i] * c_exps[i] + 4.0 * g_yy_zz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_y_z[i] = -2.0 * g_0_zz_yz_z[i] * c_exps[i] + 4.0 * g_yy_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (915-918)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_0_zz_zz_x, g_0_zz_zz_y, g_0_zz_zz_z, g_y_0_z_0_y_zz_z_x, g_y_0_z_0_y_zz_z_y, g_y_0_z_0_y_zz_z_z, g_yy_zz_0_x, g_yy_zz_0_y, g_yy_zz_0_z, g_yy_zz_zz_x, g_yy_zz_zz_y, g_yy_zz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_zz_z_x[i] = g_0_zz_0_x[i] - 2.0 * g_0_zz_zz_x[i] * c_exps[i] - 2.0 * g_yy_zz_0_x[i] * a_exp + 4.0 * g_yy_zz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_z_y[i] = g_0_zz_0_y[i] - 2.0 * g_0_zz_zz_y[i] * c_exps[i] - 2.0 * g_yy_zz_0_y[i] * a_exp + 4.0 * g_yy_zz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_z_z[i] = g_0_zz_0_z[i] - 2.0 * g_0_zz_zz_z[i] * c_exps[i] - 2.0 * g_yy_zz_0_z[i] * a_exp + 4.0 * g_yy_zz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (918-921)

    #pragma omp simd aligned(g_y_0_z_0_z_xx_x_x, g_y_0_z_0_z_xx_x_y, g_y_0_z_0_z_xx_x_z, g_yz_xx_xz_x, g_yz_xx_xz_y, g_yz_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xx_x_x[i] = 4.0 * g_yz_xx_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_x_y[i] = 4.0 * g_yz_xx_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_x_z[i] = 4.0 * g_yz_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (921-924)

    #pragma omp simd aligned(g_y_0_z_0_z_xx_y_x, g_y_0_z_0_z_xx_y_y, g_y_0_z_0_z_xx_y_z, g_yz_xx_yz_x, g_yz_xx_yz_y, g_yz_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xx_y_x[i] = 4.0 * g_yz_xx_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_y_y[i] = 4.0 * g_yz_xx_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_y_z[i] = 4.0 * g_yz_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (924-927)

    #pragma omp simd aligned(g_y_0_z_0_z_xx_z_x, g_y_0_z_0_z_xx_z_y, g_y_0_z_0_z_xx_z_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z, g_yz_xx_zz_x, g_yz_xx_zz_y, g_yz_xx_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xx_z_x[i] = -2.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_yz_xx_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_z_y[i] = -2.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_yz_xx_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_z_z[i] = -2.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_yz_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (927-930)

    #pragma omp simd aligned(g_y_0_z_0_z_xy_x_x, g_y_0_z_0_z_xy_x_y, g_y_0_z_0_z_xy_x_z, g_yz_xy_xz_x, g_yz_xy_xz_y, g_yz_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xy_x_x[i] = 4.0 * g_yz_xy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_x_y[i] = 4.0 * g_yz_xy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_x_z[i] = 4.0 * g_yz_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (930-933)

    #pragma omp simd aligned(g_y_0_z_0_z_xy_y_x, g_y_0_z_0_z_xy_y_y, g_y_0_z_0_z_xy_y_z, g_yz_xy_yz_x, g_yz_xy_yz_y, g_yz_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xy_y_x[i] = 4.0 * g_yz_xy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_y_y[i] = 4.0 * g_yz_xy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_y_z[i] = 4.0 * g_yz_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (933-936)

    #pragma omp simd aligned(g_y_0_z_0_z_xy_z_x, g_y_0_z_0_z_xy_z_y, g_y_0_z_0_z_xy_z_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_yz_xy_zz_x, g_yz_xy_zz_y, g_yz_xy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xy_z_x[i] = -2.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_yz_xy_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_z_y[i] = -2.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_yz_xy_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_z_z[i] = -2.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_yz_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (936-939)

    #pragma omp simd aligned(g_y_0_z_0_z_xz_x_x, g_y_0_z_0_z_xz_x_y, g_y_0_z_0_z_xz_x_z, g_yz_xz_xz_x, g_yz_xz_xz_y, g_yz_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xz_x_x[i] = 4.0 * g_yz_xz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_x_y[i] = 4.0 * g_yz_xz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_x_z[i] = 4.0 * g_yz_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (939-942)

    #pragma omp simd aligned(g_y_0_z_0_z_xz_y_x, g_y_0_z_0_z_xz_y_y, g_y_0_z_0_z_xz_y_z, g_yz_xz_yz_x, g_yz_xz_yz_y, g_yz_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xz_y_x[i] = 4.0 * g_yz_xz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_y_y[i] = 4.0 * g_yz_xz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_y_z[i] = 4.0 * g_yz_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (942-945)

    #pragma omp simd aligned(g_y_0_z_0_z_xz_z_x, g_y_0_z_0_z_xz_z_y, g_y_0_z_0_z_xz_z_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_yz_xz_zz_x, g_yz_xz_zz_y, g_yz_xz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xz_z_x[i] = -2.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_yz_xz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_z_y[i] = -2.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_yz_xz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_z_z[i] = -2.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_yz_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (945-948)

    #pragma omp simd aligned(g_y_0_z_0_z_yy_x_x, g_y_0_z_0_z_yy_x_y, g_y_0_z_0_z_yy_x_z, g_yz_yy_xz_x, g_yz_yy_xz_y, g_yz_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_yy_x_x[i] = 4.0 * g_yz_yy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_x_y[i] = 4.0 * g_yz_yy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_x_z[i] = 4.0 * g_yz_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (948-951)

    #pragma omp simd aligned(g_y_0_z_0_z_yy_y_x, g_y_0_z_0_z_yy_y_y, g_y_0_z_0_z_yy_y_z, g_yz_yy_yz_x, g_yz_yy_yz_y, g_yz_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_yy_y_x[i] = 4.0 * g_yz_yy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_y_y[i] = 4.0 * g_yz_yy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_y_z[i] = 4.0 * g_yz_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (951-954)

    #pragma omp simd aligned(g_y_0_z_0_z_yy_z_x, g_y_0_z_0_z_yy_z_y, g_y_0_z_0_z_yy_z_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z, g_yz_yy_zz_x, g_yz_yy_zz_y, g_yz_yy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_yy_z_x[i] = -2.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_yz_yy_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_z_y[i] = -2.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_yz_yy_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_z_z[i] = -2.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_yz_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (954-957)

    #pragma omp simd aligned(g_y_0_z_0_z_yz_x_x, g_y_0_z_0_z_yz_x_y, g_y_0_z_0_z_yz_x_z, g_yz_yz_xz_x, g_yz_yz_xz_y, g_yz_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_yz_x_x[i] = 4.0 * g_yz_yz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_x_y[i] = 4.0 * g_yz_yz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_x_z[i] = 4.0 * g_yz_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (957-960)

    #pragma omp simd aligned(g_y_0_z_0_z_yz_y_x, g_y_0_z_0_z_yz_y_y, g_y_0_z_0_z_yz_y_z, g_yz_yz_yz_x, g_yz_yz_yz_y, g_yz_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_yz_y_x[i] = 4.0 * g_yz_yz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_y_y[i] = 4.0 * g_yz_yz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_y_z[i] = 4.0 * g_yz_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (960-963)

    #pragma omp simd aligned(g_y_0_z_0_z_yz_z_x, g_y_0_z_0_z_yz_z_y, g_y_0_z_0_z_yz_z_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_yz_yz_zz_x, g_yz_yz_zz_y, g_yz_yz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_yz_z_x[i] = -2.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_yz_yz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_z_y[i] = -2.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_yz_yz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_z_z[i] = -2.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_yz_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (963-966)

    #pragma omp simd aligned(g_y_0_z_0_z_zz_x_x, g_y_0_z_0_z_zz_x_y, g_y_0_z_0_z_zz_x_z, g_yz_zz_xz_x, g_yz_zz_xz_y, g_yz_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_zz_x_x[i] = 4.0 * g_yz_zz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_x_y[i] = 4.0 * g_yz_zz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_x_z[i] = 4.0 * g_yz_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (966-969)

    #pragma omp simd aligned(g_y_0_z_0_z_zz_y_x, g_y_0_z_0_z_zz_y_y, g_y_0_z_0_z_zz_y_z, g_yz_zz_yz_x, g_yz_zz_yz_y, g_yz_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_zz_y_x[i] = 4.0 * g_yz_zz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_y_y[i] = 4.0 * g_yz_zz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_y_z[i] = 4.0 * g_yz_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (969-972)

    #pragma omp simd aligned(g_y_0_z_0_z_zz_z_x, g_y_0_z_0_z_zz_z_y, g_y_0_z_0_z_zz_z_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z, g_yz_zz_zz_x, g_yz_zz_zz_y, g_yz_zz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_zz_z_x[i] = -2.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_yz_zz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_z_y[i] = -2.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_yz_zz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_z_z[i] = -2.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_yz_zz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (972-975)

    #pragma omp simd aligned(g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_xz_xx_xx_x, g_xz_xx_xx_y, g_xz_xx_xx_z, g_z_0_x_0_x_xx_x_x, g_z_0_x_0_x_xx_x_y, g_z_0_x_0_x_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xx_x_x[i] = -2.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xz_xx_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_x_y[i] = -2.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xz_xx_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_x_z[i] = -2.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xz_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (975-978)

    #pragma omp simd aligned(g_xz_xx_xy_x, g_xz_xx_xy_y, g_xz_xx_xy_z, g_z_0_x_0_x_xx_y_x, g_z_0_x_0_x_xx_y_y, g_z_0_x_0_x_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xx_y_x[i] = 4.0 * g_xz_xx_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_y_y[i] = 4.0 * g_xz_xx_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_y_z[i] = 4.0 * g_xz_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (978-981)

    #pragma omp simd aligned(g_xz_xx_xz_x, g_xz_xx_xz_y, g_xz_xx_xz_z, g_z_0_x_0_x_xx_z_x, g_z_0_x_0_x_xx_z_y, g_z_0_x_0_x_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xx_z_x[i] = 4.0 * g_xz_xx_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_z_y[i] = 4.0 * g_xz_xx_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_z_z[i] = 4.0 * g_xz_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (981-984)

    #pragma omp simd aligned(g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_xz_xy_xx_x, g_xz_xy_xx_y, g_xz_xy_xx_z, g_z_0_x_0_x_xy_x_x, g_z_0_x_0_x_xy_x_y, g_z_0_x_0_x_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xy_x_x[i] = -2.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xz_xy_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_x_y[i] = -2.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xz_xy_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_x_z[i] = -2.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xz_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (984-987)

    #pragma omp simd aligned(g_xz_xy_xy_x, g_xz_xy_xy_y, g_xz_xy_xy_z, g_z_0_x_0_x_xy_y_x, g_z_0_x_0_x_xy_y_y, g_z_0_x_0_x_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xy_y_x[i] = 4.0 * g_xz_xy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_y_y[i] = 4.0 * g_xz_xy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_y_z[i] = 4.0 * g_xz_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (987-990)

    #pragma omp simd aligned(g_xz_xy_xz_x, g_xz_xy_xz_y, g_xz_xy_xz_z, g_z_0_x_0_x_xy_z_x, g_z_0_x_0_x_xy_z_y, g_z_0_x_0_x_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xy_z_x[i] = 4.0 * g_xz_xy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_z_y[i] = 4.0 * g_xz_xy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_z_z[i] = 4.0 * g_xz_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (990-993)

    #pragma omp simd aligned(g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_xz_xz_xx_x, g_xz_xz_xx_y, g_xz_xz_xx_z, g_z_0_x_0_x_xz_x_x, g_z_0_x_0_x_xz_x_y, g_z_0_x_0_x_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xz_x_x[i] = -2.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xz_xz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_x_y[i] = -2.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xz_xz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_x_z[i] = -2.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xz_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (993-996)

    #pragma omp simd aligned(g_xz_xz_xy_x, g_xz_xz_xy_y, g_xz_xz_xy_z, g_z_0_x_0_x_xz_y_x, g_z_0_x_0_x_xz_y_y, g_z_0_x_0_x_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xz_y_x[i] = 4.0 * g_xz_xz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_y_y[i] = 4.0 * g_xz_xz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_y_z[i] = 4.0 * g_xz_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (996-999)

    #pragma omp simd aligned(g_xz_xz_xz_x, g_xz_xz_xz_y, g_xz_xz_xz_z, g_z_0_x_0_x_xz_z_x, g_z_0_x_0_x_xz_z_y, g_z_0_x_0_x_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xz_z_x[i] = 4.0 * g_xz_xz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_z_y[i] = 4.0 * g_xz_xz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_z_z[i] = 4.0 * g_xz_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (999-1002)

    #pragma omp simd aligned(g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_xz_yy_xx_x, g_xz_yy_xx_y, g_xz_yy_xx_z, g_z_0_x_0_x_yy_x_x, g_z_0_x_0_x_yy_x_y, g_z_0_x_0_x_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_yy_x_x[i] = -2.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xz_yy_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_x_y[i] = -2.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xz_yy_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_x_z[i] = -2.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xz_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1002-1005)

    #pragma omp simd aligned(g_xz_yy_xy_x, g_xz_yy_xy_y, g_xz_yy_xy_z, g_z_0_x_0_x_yy_y_x, g_z_0_x_0_x_yy_y_y, g_z_0_x_0_x_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_yy_y_x[i] = 4.0 * g_xz_yy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_y_y[i] = 4.0 * g_xz_yy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_y_z[i] = 4.0 * g_xz_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1005-1008)

    #pragma omp simd aligned(g_xz_yy_xz_x, g_xz_yy_xz_y, g_xz_yy_xz_z, g_z_0_x_0_x_yy_z_x, g_z_0_x_0_x_yy_z_y, g_z_0_x_0_x_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_yy_z_x[i] = 4.0 * g_xz_yy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_z_y[i] = 4.0 * g_xz_yy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_z_z[i] = 4.0 * g_xz_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1008-1011)

    #pragma omp simd aligned(g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_xz_yz_xx_x, g_xz_yz_xx_y, g_xz_yz_xx_z, g_z_0_x_0_x_yz_x_x, g_z_0_x_0_x_yz_x_y, g_z_0_x_0_x_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_yz_x_x[i] = -2.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xz_yz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_x_y[i] = -2.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xz_yz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_x_z[i] = -2.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xz_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1011-1014)

    #pragma omp simd aligned(g_xz_yz_xy_x, g_xz_yz_xy_y, g_xz_yz_xy_z, g_z_0_x_0_x_yz_y_x, g_z_0_x_0_x_yz_y_y, g_z_0_x_0_x_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_yz_y_x[i] = 4.0 * g_xz_yz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_y_y[i] = 4.0 * g_xz_yz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_y_z[i] = 4.0 * g_xz_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1014-1017)

    #pragma omp simd aligned(g_xz_yz_xz_x, g_xz_yz_xz_y, g_xz_yz_xz_z, g_z_0_x_0_x_yz_z_x, g_z_0_x_0_x_yz_z_y, g_z_0_x_0_x_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_yz_z_x[i] = 4.0 * g_xz_yz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_z_y[i] = 4.0 * g_xz_yz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_z_z[i] = 4.0 * g_xz_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1017-1020)

    #pragma omp simd aligned(g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_xz_zz_xx_x, g_xz_zz_xx_y, g_xz_zz_xx_z, g_z_0_x_0_x_zz_x_x, g_z_0_x_0_x_zz_x_y, g_z_0_x_0_x_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_zz_x_x[i] = -2.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xz_zz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_x_y[i] = -2.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xz_zz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_x_z[i] = -2.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xz_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1020-1023)

    #pragma omp simd aligned(g_xz_zz_xy_x, g_xz_zz_xy_y, g_xz_zz_xy_z, g_z_0_x_0_x_zz_y_x, g_z_0_x_0_x_zz_y_y, g_z_0_x_0_x_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_zz_y_x[i] = 4.0 * g_xz_zz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_y_y[i] = 4.0 * g_xz_zz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_y_z[i] = 4.0 * g_xz_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1023-1026)

    #pragma omp simd aligned(g_xz_zz_xz_x, g_xz_zz_xz_y, g_xz_zz_xz_z, g_z_0_x_0_x_zz_z_x, g_z_0_x_0_x_zz_z_y, g_z_0_x_0_x_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_zz_z_x[i] = 4.0 * g_xz_zz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_z_y[i] = 4.0 * g_xz_zz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_z_z[i] = 4.0 * g_xz_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1026-1029)

    #pragma omp simd aligned(g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z, g_yz_xx_xx_x, g_yz_xx_xx_y, g_yz_xx_xx_z, g_z_0_x_0_y_xx_x_x, g_z_0_x_0_y_xx_x_y, g_z_0_x_0_y_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xx_x_x[i] = -2.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_yz_xx_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_x_y[i] = -2.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_yz_xx_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_x_z[i] = -2.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_yz_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1029-1032)

    #pragma omp simd aligned(g_yz_xx_xy_x, g_yz_xx_xy_y, g_yz_xx_xy_z, g_z_0_x_0_y_xx_y_x, g_z_0_x_0_y_xx_y_y, g_z_0_x_0_y_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xx_y_x[i] = 4.0 * g_yz_xx_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_y_y[i] = 4.0 * g_yz_xx_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_y_z[i] = 4.0 * g_yz_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1032-1035)

    #pragma omp simd aligned(g_yz_xx_xz_x, g_yz_xx_xz_y, g_yz_xx_xz_z, g_z_0_x_0_y_xx_z_x, g_z_0_x_0_y_xx_z_y, g_z_0_x_0_y_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xx_z_x[i] = 4.0 * g_yz_xx_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_z_y[i] = 4.0 * g_yz_xx_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_z_z[i] = 4.0 * g_yz_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1035-1038)

    #pragma omp simd aligned(g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_yz_xy_xx_x, g_yz_xy_xx_y, g_yz_xy_xx_z, g_z_0_x_0_y_xy_x_x, g_z_0_x_0_y_xy_x_y, g_z_0_x_0_y_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xy_x_x[i] = -2.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_yz_xy_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_x_y[i] = -2.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_yz_xy_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_x_z[i] = -2.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_yz_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1038-1041)

    #pragma omp simd aligned(g_yz_xy_xy_x, g_yz_xy_xy_y, g_yz_xy_xy_z, g_z_0_x_0_y_xy_y_x, g_z_0_x_0_y_xy_y_y, g_z_0_x_0_y_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xy_y_x[i] = 4.0 * g_yz_xy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_y_y[i] = 4.0 * g_yz_xy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_y_z[i] = 4.0 * g_yz_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1041-1044)

    #pragma omp simd aligned(g_yz_xy_xz_x, g_yz_xy_xz_y, g_yz_xy_xz_z, g_z_0_x_0_y_xy_z_x, g_z_0_x_0_y_xy_z_y, g_z_0_x_0_y_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xy_z_x[i] = 4.0 * g_yz_xy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_z_y[i] = 4.0 * g_yz_xy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_z_z[i] = 4.0 * g_yz_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1044-1047)

    #pragma omp simd aligned(g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_yz_xz_xx_x, g_yz_xz_xx_y, g_yz_xz_xx_z, g_z_0_x_0_y_xz_x_x, g_z_0_x_0_y_xz_x_y, g_z_0_x_0_y_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xz_x_x[i] = -2.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_yz_xz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_x_y[i] = -2.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_yz_xz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_x_z[i] = -2.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_yz_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1047-1050)

    #pragma omp simd aligned(g_yz_xz_xy_x, g_yz_xz_xy_y, g_yz_xz_xy_z, g_z_0_x_0_y_xz_y_x, g_z_0_x_0_y_xz_y_y, g_z_0_x_0_y_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xz_y_x[i] = 4.0 * g_yz_xz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_y_y[i] = 4.0 * g_yz_xz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_y_z[i] = 4.0 * g_yz_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1050-1053)

    #pragma omp simd aligned(g_yz_xz_xz_x, g_yz_xz_xz_y, g_yz_xz_xz_z, g_z_0_x_0_y_xz_z_x, g_z_0_x_0_y_xz_z_y, g_z_0_x_0_y_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xz_z_x[i] = 4.0 * g_yz_xz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_z_y[i] = 4.0 * g_yz_xz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_z_z[i] = 4.0 * g_yz_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1053-1056)

    #pragma omp simd aligned(g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z, g_yz_yy_xx_x, g_yz_yy_xx_y, g_yz_yy_xx_z, g_z_0_x_0_y_yy_x_x, g_z_0_x_0_y_yy_x_y, g_z_0_x_0_y_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_yy_x_x[i] = -2.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_yz_yy_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_x_y[i] = -2.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_yz_yy_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_x_z[i] = -2.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_yz_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1056-1059)

    #pragma omp simd aligned(g_yz_yy_xy_x, g_yz_yy_xy_y, g_yz_yy_xy_z, g_z_0_x_0_y_yy_y_x, g_z_0_x_0_y_yy_y_y, g_z_0_x_0_y_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_yy_y_x[i] = 4.0 * g_yz_yy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_y_y[i] = 4.0 * g_yz_yy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_y_z[i] = 4.0 * g_yz_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1059-1062)

    #pragma omp simd aligned(g_yz_yy_xz_x, g_yz_yy_xz_y, g_yz_yy_xz_z, g_z_0_x_0_y_yy_z_x, g_z_0_x_0_y_yy_z_y, g_z_0_x_0_y_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_yy_z_x[i] = 4.0 * g_yz_yy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_z_y[i] = 4.0 * g_yz_yy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_z_z[i] = 4.0 * g_yz_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1062-1065)

    #pragma omp simd aligned(g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_yz_yz_xx_x, g_yz_yz_xx_y, g_yz_yz_xx_z, g_z_0_x_0_y_yz_x_x, g_z_0_x_0_y_yz_x_y, g_z_0_x_0_y_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_yz_x_x[i] = -2.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_yz_yz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_x_y[i] = -2.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_yz_yz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_x_z[i] = -2.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_yz_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1065-1068)

    #pragma omp simd aligned(g_yz_yz_xy_x, g_yz_yz_xy_y, g_yz_yz_xy_z, g_z_0_x_0_y_yz_y_x, g_z_0_x_0_y_yz_y_y, g_z_0_x_0_y_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_yz_y_x[i] = 4.0 * g_yz_yz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_y_y[i] = 4.0 * g_yz_yz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_y_z[i] = 4.0 * g_yz_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1068-1071)

    #pragma omp simd aligned(g_yz_yz_xz_x, g_yz_yz_xz_y, g_yz_yz_xz_z, g_z_0_x_0_y_yz_z_x, g_z_0_x_0_y_yz_z_y, g_z_0_x_0_y_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_yz_z_x[i] = 4.0 * g_yz_yz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_z_y[i] = 4.0 * g_yz_yz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_z_z[i] = 4.0 * g_yz_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1071-1074)

    #pragma omp simd aligned(g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z, g_yz_zz_xx_x, g_yz_zz_xx_y, g_yz_zz_xx_z, g_z_0_x_0_y_zz_x_x, g_z_0_x_0_y_zz_x_y, g_z_0_x_0_y_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_zz_x_x[i] = -2.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_yz_zz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_x_y[i] = -2.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_yz_zz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_x_z[i] = -2.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_yz_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1074-1077)

    #pragma omp simd aligned(g_yz_zz_xy_x, g_yz_zz_xy_y, g_yz_zz_xy_z, g_z_0_x_0_y_zz_y_x, g_z_0_x_0_y_zz_y_y, g_z_0_x_0_y_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_zz_y_x[i] = 4.0 * g_yz_zz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_y_y[i] = 4.0 * g_yz_zz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_y_z[i] = 4.0 * g_yz_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1077-1080)

    #pragma omp simd aligned(g_yz_zz_xz_x, g_yz_zz_xz_y, g_yz_zz_xz_z, g_z_0_x_0_y_zz_z_x, g_z_0_x_0_y_zz_z_y, g_z_0_x_0_y_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_zz_z_x[i] = 4.0 * g_yz_zz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_z_y[i] = 4.0 * g_yz_zz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_z_z[i] = 4.0 * g_yz_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1080-1083)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_0_xx_xx_x, g_0_xx_xx_y, g_0_xx_xx_z, g_z_0_x_0_z_xx_x_x, g_z_0_x_0_z_xx_x_y, g_z_0_x_0_z_xx_x_z, g_zz_xx_0_x, g_zz_xx_0_y, g_zz_xx_0_z, g_zz_xx_xx_x, g_zz_xx_xx_y, g_zz_xx_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xx_x_x[i] = g_0_xx_0_x[i] - 2.0 * g_0_xx_xx_x[i] * c_exps[i] - 2.0 * g_zz_xx_0_x[i] * a_exp + 4.0 * g_zz_xx_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_x_y[i] = g_0_xx_0_y[i] - 2.0 * g_0_xx_xx_y[i] * c_exps[i] - 2.0 * g_zz_xx_0_y[i] * a_exp + 4.0 * g_zz_xx_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_x_z[i] = g_0_xx_0_z[i] - 2.0 * g_0_xx_xx_z[i] * c_exps[i] - 2.0 * g_zz_xx_0_z[i] * a_exp + 4.0 * g_zz_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1083-1086)

    #pragma omp simd aligned(g_0_xx_xy_x, g_0_xx_xy_y, g_0_xx_xy_z, g_z_0_x_0_z_xx_y_x, g_z_0_x_0_z_xx_y_y, g_z_0_x_0_z_xx_y_z, g_zz_xx_xy_x, g_zz_xx_xy_y, g_zz_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xx_y_x[i] = -2.0 * g_0_xx_xy_x[i] * c_exps[i] + 4.0 * g_zz_xx_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_y_y[i] = -2.0 * g_0_xx_xy_y[i] * c_exps[i] + 4.0 * g_zz_xx_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_y_z[i] = -2.0 * g_0_xx_xy_z[i] * c_exps[i] + 4.0 * g_zz_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1086-1089)

    #pragma omp simd aligned(g_0_xx_xz_x, g_0_xx_xz_y, g_0_xx_xz_z, g_z_0_x_0_z_xx_z_x, g_z_0_x_0_z_xx_z_y, g_z_0_x_0_z_xx_z_z, g_zz_xx_xz_x, g_zz_xx_xz_y, g_zz_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xx_z_x[i] = -2.0 * g_0_xx_xz_x[i] * c_exps[i] + 4.0 * g_zz_xx_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_z_y[i] = -2.0 * g_0_xx_xz_y[i] * c_exps[i] + 4.0 * g_zz_xx_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_z_z[i] = -2.0 * g_0_xx_xz_z[i] * c_exps[i] + 4.0 * g_zz_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1089-1092)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_0_xy_xx_x, g_0_xy_xx_y, g_0_xy_xx_z, g_z_0_x_0_z_xy_x_x, g_z_0_x_0_z_xy_x_y, g_z_0_x_0_z_xy_x_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z, g_zz_xy_xx_x, g_zz_xy_xx_y, g_zz_xy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xy_x_x[i] = g_0_xy_0_x[i] - 2.0 * g_0_xy_xx_x[i] * c_exps[i] - 2.0 * g_zz_xy_0_x[i] * a_exp + 4.0 * g_zz_xy_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_x_y[i] = g_0_xy_0_y[i] - 2.0 * g_0_xy_xx_y[i] * c_exps[i] - 2.0 * g_zz_xy_0_y[i] * a_exp + 4.0 * g_zz_xy_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_x_z[i] = g_0_xy_0_z[i] - 2.0 * g_0_xy_xx_z[i] * c_exps[i] - 2.0 * g_zz_xy_0_z[i] * a_exp + 4.0 * g_zz_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1092-1095)

    #pragma omp simd aligned(g_0_xy_xy_x, g_0_xy_xy_y, g_0_xy_xy_z, g_z_0_x_0_z_xy_y_x, g_z_0_x_0_z_xy_y_y, g_z_0_x_0_z_xy_y_z, g_zz_xy_xy_x, g_zz_xy_xy_y, g_zz_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xy_y_x[i] = -2.0 * g_0_xy_xy_x[i] * c_exps[i] + 4.0 * g_zz_xy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_y_y[i] = -2.0 * g_0_xy_xy_y[i] * c_exps[i] + 4.0 * g_zz_xy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_y_z[i] = -2.0 * g_0_xy_xy_z[i] * c_exps[i] + 4.0 * g_zz_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1095-1098)

    #pragma omp simd aligned(g_0_xy_xz_x, g_0_xy_xz_y, g_0_xy_xz_z, g_z_0_x_0_z_xy_z_x, g_z_0_x_0_z_xy_z_y, g_z_0_x_0_z_xy_z_z, g_zz_xy_xz_x, g_zz_xy_xz_y, g_zz_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xy_z_x[i] = -2.0 * g_0_xy_xz_x[i] * c_exps[i] + 4.0 * g_zz_xy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_z_y[i] = -2.0 * g_0_xy_xz_y[i] * c_exps[i] + 4.0 * g_zz_xy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_z_z[i] = -2.0 * g_0_xy_xz_z[i] * c_exps[i] + 4.0 * g_zz_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1098-1101)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_0_xz_xx_x, g_0_xz_xx_y, g_0_xz_xx_z, g_z_0_x_0_z_xz_x_x, g_z_0_x_0_z_xz_x_y, g_z_0_x_0_z_xz_x_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z, g_zz_xz_xx_x, g_zz_xz_xx_y, g_zz_xz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xz_x_x[i] = g_0_xz_0_x[i] - 2.0 * g_0_xz_xx_x[i] * c_exps[i] - 2.0 * g_zz_xz_0_x[i] * a_exp + 4.0 * g_zz_xz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_x_y[i] = g_0_xz_0_y[i] - 2.0 * g_0_xz_xx_y[i] * c_exps[i] - 2.0 * g_zz_xz_0_y[i] * a_exp + 4.0 * g_zz_xz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_x_z[i] = g_0_xz_0_z[i] - 2.0 * g_0_xz_xx_z[i] * c_exps[i] - 2.0 * g_zz_xz_0_z[i] * a_exp + 4.0 * g_zz_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1101-1104)

    #pragma omp simd aligned(g_0_xz_xy_x, g_0_xz_xy_y, g_0_xz_xy_z, g_z_0_x_0_z_xz_y_x, g_z_0_x_0_z_xz_y_y, g_z_0_x_0_z_xz_y_z, g_zz_xz_xy_x, g_zz_xz_xy_y, g_zz_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xz_y_x[i] = -2.0 * g_0_xz_xy_x[i] * c_exps[i] + 4.0 * g_zz_xz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_y_y[i] = -2.0 * g_0_xz_xy_y[i] * c_exps[i] + 4.0 * g_zz_xz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_y_z[i] = -2.0 * g_0_xz_xy_z[i] * c_exps[i] + 4.0 * g_zz_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1104-1107)

    #pragma omp simd aligned(g_0_xz_xz_x, g_0_xz_xz_y, g_0_xz_xz_z, g_z_0_x_0_z_xz_z_x, g_z_0_x_0_z_xz_z_y, g_z_0_x_0_z_xz_z_z, g_zz_xz_xz_x, g_zz_xz_xz_y, g_zz_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xz_z_x[i] = -2.0 * g_0_xz_xz_x[i] * c_exps[i] + 4.0 * g_zz_xz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_z_y[i] = -2.0 * g_0_xz_xz_y[i] * c_exps[i] + 4.0 * g_zz_xz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_z_z[i] = -2.0 * g_0_xz_xz_z[i] * c_exps[i] + 4.0 * g_zz_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1107-1110)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_0_yy_xx_x, g_0_yy_xx_y, g_0_yy_xx_z, g_z_0_x_0_z_yy_x_x, g_z_0_x_0_z_yy_x_y, g_z_0_x_0_z_yy_x_z, g_zz_yy_0_x, g_zz_yy_0_y, g_zz_yy_0_z, g_zz_yy_xx_x, g_zz_yy_xx_y, g_zz_yy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_yy_x_x[i] = g_0_yy_0_x[i] - 2.0 * g_0_yy_xx_x[i] * c_exps[i] - 2.0 * g_zz_yy_0_x[i] * a_exp + 4.0 * g_zz_yy_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_x_y[i] = g_0_yy_0_y[i] - 2.0 * g_0_yy_xx_y[i] * c_exps[i] - 2.0 * g_zz_yy_0_y[i] * a_exp + 4.0 * g_zz_yy_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_x_z[i] = g_0_yy_0_z[i] - 2.0 * g_0_yy_xx_z[i] * c_exps[i] - 2.0 * g_zz_yy_0_z[i] * a_exp + 4.0 * g_zz_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1110-1113)

    #pragma omp simd aligned(g_0_yy_xy_x, g_0_yy_xy_y, g_0_yy_xy_z, g_z_0_x_0_z_yy_y_x, g_z_0_x_0_z_yy_y_y, g_z_0_x_0_z_yy_y_z, g_zz_yy_xy_x, g_zz_yy_xy_y, g_zz_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_yy_y_x[i] = -2.0 * g_0_yy_xy_x[i] * c_exps[i] + 4.0 * g_zz_yy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_y_y[i] = -2.0 * g_0_yy_xy_y[i] * c_exps[i] + 4.0 * g_zz_yy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_y_z[i] = -2.0 * g_0_yy_xy_z[i] * c_exps[i] + 4.0 * g_zz_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1113-1116)

    #pragma omp simd aligned(g_0_yy_xz_x, g_0_yy_xz_y, g_0_yy_xz_z, g_z_0_x_0_z_yy_z_x, g_z_0_x_0_z_yy_z_y, g_z_0_x_0_z_yy_z_z, g_zz_yy_xz_x, g_zz_yy_xz_y, g_zz_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_yy_z_x[i] = -2.0 * g_0_yy_xz_x[i] * c_exps[i] + 4.0 * g_zz_yy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_z_y[i] = -2.0 * g_0_yy_xz_y[i] * c_exps[i] + 4.0 * g_zz_yy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_z_z[i] = -2.0 * g_0_yy_xz_z[i] * c_exps[i] + 4.0 * g_zz_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1116-1119)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_0_yz_xx_x, g_0_yz_xx_y, g_0_yz_xx_z, g_z_0_x_0_z_yz_x_x, g_z_0_x_0_z_yz_x_y, g_z_0_x_0_z_yz_x_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z, g_zz_yz_xx_x, g_zz_yz_xx_y, g_zz_yz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_yz_x_x[i] = g_0_yz_0_x[i] - 2.0 * g_0_yz_xx_x[i] * c_exps[i] - 2.0 * g_zz_yz_0_x[i] * a_exp + 4.0 * g_zz_yz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_x_y[i] = g_0_yz_0_y[i] - 2.0 * g_0_yz_xx_y[i] * c_exps[i] - 2.0 * g_zz_yz_0_y[i] * a_exp + 4.0 * g_zz_yz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_x_z[i] = g_0_yz_0_z[i] - 2.0 * g_0_yz_xx_z[i] * c_exps[i] - 2.0 * g_zz_yz_0_z[i] * a_exp + 4.0 * g_zz_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1119-1122)

    #pragma omp simd aligned(g_0_yz_xy_x, g_0_yz_xy_y, g_0_yz_xy_z, g_z_0_x_0_z_yz_y_x, g_z_0_x_0_z_yz_y_y, g_z_0_x_0_z_yz_y_z, g_zz_yz_xy_x, g_zz_yz_xy_y, g_zz_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_yz_y_x[i] = -2.0 * g_0_yz_xy_x[i] * c_exps[i] + 4.0 * g_zz_yz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_y_y[i] = -2.0 * g_0_yz_xy_y[i] * c_exps[i] + 4.0 * g_zz_yz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_y_z[i] = -2.0 * g_0_yz_xy_z[i] * c_exps[i] + 4.0 * g_zz_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1122-1125)

    #pragma omp simd aligned(g_0_yz_xz_x, g_0_yz_xz_y, g_0_yz_xz_z, g_z_0_x_0_z_yz_z_x, g_z_0_x_0_z_yz_z_y, g_z_0_x_0_z_yz_z_z, g_zz_yz_xz_x, g_zz_yz_xz_y, g_zz_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_yz_z_x[i] = -2.0 * g_0_yz_xz_x[i] * c_exps[i] + 4.0 * g_zz_yz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_z_y[i] = -2.0 * g_0_yz_xz_y[i] * c_exps[i] + 4.0 * g_zz_yz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_z_z[i] = -2.0 * g_0_yz_xz_z[i] * c_exps[i] + 4.0 * g_zz_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1125-1128)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_0_zz_xx_x, g_0_zz_xx_y, g_0_zz_xx_z, g_z_0_x_0_z_zz_x_x, g_z_0_x_0_z_zz_x_y, g_z_0_x_0_z_zz_x_z, g_zz_zz_0_x, g_zz_zz_0_y, g_zz_zz_0_z, g_zz_zz_xx_x, g_zz_zz_xx_y, g_zz_zz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_zz_x_x[i] = g_0_zz_0_x[i] - 2.0 * g_0_zz_xx_x[i] * c_exps[i] - 2.0 * g_zz_zz_0_x[i] * a_exp + 4.0 * g_zz_zz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_x_y[i] = g_0_zz_0_y[i] - 2.0 * g_0_zz_xx_y[i] * c_exps[i] - 2.0 * g_zz_zz_0_y[i] * a_exp + 4.0 * g_zz_zz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_x_z[i] = g_0_zz_0_z[i] - 2.0 * g_0_zz_xx_z[i] * c_exps[i] - 2.0 * g_zz_zz_0_z[i] * a_exp + 4.0 * g_zz_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1128-1131)

    #pragma omp simd aligned(g_0_zz_xy_x, g_0_zz_xy_y, g_0_zz_xy_z, g_z_0_x_0_z_zz_y_x, g_z_0_x_0_z_zz_y_y, g_z_0_x_0_z_zz_y_z, g_zz_zz_xy_x, g_zz_zz_xy_y, g_zz_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_zz_y_x[i] = -2.0 * g_0_zz_xy_x[i] * c_exps[i] + 4.0 * g_zz_zz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_y_y[i] = -2.0 * g_0_zz_xy_y[i] * c_exps[i] + 4.0 * g_zz_zz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_y_z[i] = -2.0 * g_0_zz_xy_z[i] * c_exps[i] + 4.0 * g_zz_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1131-1134)

    #pragma omp simd aligned(g_0_zz_xz_x, g_0_zz_xz_y, g_0_zz_xz_z, g_z_0_x_0_z_zz_z_x, g_z_0_x_0_z_zz_z_y, g_z_0_x_0_z_zz_z_z, g_zz_zz_xz_x, g_zz_zz_xz_y, g_zz_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_zz_z_x[i] = -2.0 * g_0_zz_xz_x[i] * c_exps[i] + 4.0 * g_zz_zz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_z_y[i] = -2.0 * g_0_zz_xz_y[i] * c_exps[i] + 4.0 * g_zz_zz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_z_z[i] = -2.0 * g_0_zz_xz_z[i] * c_exps[i] + 4.0 * g_zz_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1134-1137)

    #pragma omp simd aligned(g_xz_xx_xy_x, g_xz_xx_xy_y, g_xz_xx_xy_z, g_z_0_y_0_x_xx_x_x, g_z_0_y_0_x_xx_x_y, g_z_0_y_0_x_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xx_x_x[i] = 4.0 * g_xz_xx_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_x_y[i] = 4.0 * g_xz_xx_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_x_z[i] = 4.0 * g_xz_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1137-1140)

    #pragma omp simd aligned(g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_xz_xx_yy_x, g_xz_xx_yy_y, g_xz_xx_yy_z, g_z_0_y_0_x_xx_y_x, g_z_0_y_0_x_xx_y_y, g_z_0_y_0_x_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xx_y_x[i] = -2.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xz_xx_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_y_y[i] = -2.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xz_xx_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_y_z[i] = -2.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xz_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1140-1143)

    #pragma omp simd aligned(g_xz_xx_yz_x, g_xz_xx_yz_y, g_xz_xx_yz_z, g_z_0_y_0_x_xx_z_x, g_z_0_y_0_x_xx_z_y, g_z_0_y_0_x_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xx_z_x[i] = 4.0 * g_xz_xx_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_z_y[i] = 4.0 * g_xz_xx_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_z_z[i] = 4.0 * g_xz_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1143-1146)

    #pragma omp simd aligned(g_xz_xy_xy_x, g_xz_xy_xy_y, g_xz_xy_xy_z, g_z_0_y_0_x_xy_x_x, g_z_0_y_0_x_xy_x_y, g_z_0_y_0_x_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xy_x_x[i] = 4.0 * g_xz_xy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_x_y[i] = 4.0 * g_xz_xy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_x_z[i] = 4.0 * g_xz_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1146-1149)

    #pragma omp simd aligned(g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_xz_xy_yy_x, g_xz_xy_yy_y, g_xz_xy_yy_z, g_z_0_y_0_x_xy_y_x, g_z_0_y_0_x_xy_y_y, g_z_0_y_0_x_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xy_y_x[i] = -2.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xz_xy_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_y_y[i] = -2.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xz_xy_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_y_z[i] = -2.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xz_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1149-1152)

    #pragma omp simd aligned(g_xz_xy_yz_x, g_xz_xy_yz_y, g_xz_xy_yz_z, g_z_0_y_0_x_xy_z_x, g_z_0_y_0_x_xy_z_y, g_z_0_y_0_x_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xy_z_x[i] = 4.0 * g_xz_xy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_z_y[i] = 4.0 * g_xz_xy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_z_z[i] = 4.0 * g_xz_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1152-1155)

    #pragma omp simd aligned(g_xz_xz_xy_x, g_xz_xz_xy_y, g_xz_xz_xy_z, g_z_0_y_0_x_xz_x_x, g_z_0_y_0_x_xz_x_y, g_z_0_y_0_x_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xz_x_x[i] = 4.0 * g_xz_xz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_x_y[i] = 4.0 * g_xz_xz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_x_z[i] = 4.0 * g_xz_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1155-1158)

    #pragma omp simd aligned(g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_xz_xz_yy_x, g_xz_xz_yy_y, g_xz_xz_yy_z, g_z_0_y_0_x_xz_y_x, g_z_0_y_0_x_xz_y_y, g_z_0_y_0_x_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xz_y_x[i] = -2.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xz_xz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_y_y[i] = -2.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xz_xz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_y_z[i] = -2.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xz_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1158-1161)

    #pragma omp simd aligned(g_xz_xz_yz_x, g_xz_xz_yz_y, g_xz_xz_yz_z, g_z_0_y_0_x_xz_z_x, g_z_0_y_0_x_xz_z_y, g_z_0_y_0_x_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xz_z_x[i] = 4.0 * g_xz_xz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_z_y[i] = 4.0 * g_xz_xz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_z_z[i] = 4.0 * g_xz_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1161-1164)

    #pragma omp simd aligned(g_xz_yy_xy_x, g_xz_yy_xy_y, g_xz_yy_xy_z, g_z_0_y_0_x_yy_x_x, g_z_0_y_0_x_yy_x_y, g_z_0_y_0_x_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_yy_x_x[i] = 4.0 * g_xz_yy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_x_y[i] = 4.0 * g_xz_yy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_x_z[i] = 4.0 * g_xz_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1164-1167)

    #pragma omp simd aligned(g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_xz_yy_yy_x, g_xz_yy_yy_y, g_xz_yy_yy_z, g_z_0_y_0_x_yy_y_x, g_z_0_y_0_x_yy_y_y, g_z_0_y_0_x_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_yy_y_x[i] = -2.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xz_yy_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_y_y[i] = -2.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xz_yy_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_y_z[i] = -2.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xz_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1167-1170)

    #pragma omp simd aligned(g_xz_yy_yz_x, g_xz_yy_yz_y, g_xz_yy_yz_z, g_z_0_y_0_x_yy_z_x, g_z_0_y_0_x_yy_z_y, g_z_0_y_0_x_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_yy_z_x[i] = 4.0 * g_xz_yy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_z_y[i] = 4.0 * g_xz_yy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_z_z[i] = 4.0 * g_xz_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1170-1173)

    #pragma omp simd aligned(g_xz_yz_xy_x, g_xz_yz_xy_y, g_xz_yz_xy_z, g_z_0_y_0_x_yz_x_x, g_z_0_y_0_x_yz_x_y, g_z_0_y_0_x_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_yz_x_x[i] = 4.0 * g_xz_yz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_x_y[i] = 4.0 * g_xz_yz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_x_z[i] = 4.0 * g_xz_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1173-1176)

    #pragma omp simd aligned(g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_xz_yz_yy_x, g_xz_yz_yy_y, g_xz_yz_yy_z, g_z_0_y_0_x_yz_y_x, g_z_0_y_0_x_yz_y_y, g_z_0_y_0_x_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_yz_y_x[i] = -2.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xz_yz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_y_y[i] = -2.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xz_yz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_y_z[i] = -2.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xz_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1176-1179)

    #pragma omp simd aligned(g_xz_yz_yz_x, g_xz_yz_yz_y, g_xz_yz_yz_z, g_z_0_y_0_x_yz_z_x, g_z_0_y_0_x_yz_z_y, g_z_0_y_0_x_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_yz_z_x[i] = 4.0 * g_xz_yz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_z_y[i] = 4.0 * g_xz_yz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_z_z[i] = 4.0 * g_xz_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1179-1182)

    #pragma omp simd aligned(g_xz_zz_xy_x, g_xz_zz_xy_y, g_xz_zz_xy_z, g_z_0_y_0_x_zz_x_x, g_z_0_y_0_x_zz_x_y, g_z_0_y_0_x_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_zz_x_x[i] = 4.0 * g_xz_zz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_x_y[i] = 4.0 * g_xz_zz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_x_z[i] = 4.0 * g_xz_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1182-1185)

    #pragma omp simd aligned(g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_xz_zz_yy_x, g_xz_zz_yy_y, g_xz_zz_yy_z, g_z_0_y_0_x_zz_y_x, g_z_0_y_0_x_zz_y_y, g_z_0_y_0_x_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_zz_y_x[i] = -2.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xz_zz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_y_y[i] = -2.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xz_zz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_y_z[i] = -2.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xz_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1185-1188)

    #pragma omp simd aligned(g_xz_zz_yz_x, g_xz_zz_yz_y, g_xz_zz_yz_z, g_z_0_y_0_x_zz_z_x, g_z_0_y_0_x_zz_z_y, g_z_0_y_0_x_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_zz_z_x[i] = 4.0 * g_xz_zz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_z_y[i] = 4.0 * g_xz_zz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_z_z[i] = 4.0 * g_xz_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1188-1191)

    #pragma omp simd aligned(g_yz_xx_xy_x, g_yz_xx_xy_y, g_yz_xx_xy_z, g_z_0_y_0_y_xx_x_x, g_z_0_y_0_y_xx_x_y, g_z_0_y_0_y_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xx_x_x[i] = 4.0 * g_yz_xx_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_x_y[i] = 4.0 * g_yz_xx_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_x_z[i] = 4.0 * g_yz_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1191-1194)

    #pragma omp simd aligned(g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z, g_yz_xx_yy_x, g_yz_xx_yy_y, g_yz_xx_yy_z, g_z_0_y_0_y_xx_y_x, g_z_0_y_0_y_xx_y_y, g_z_0_y_0_y_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xx_y_x[i] = -2.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_yz_xx_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_y_y[i] = -2.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_yz_xx_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_y_z[i] = -2.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_yz_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1194-1197)

    #pragma omp simd aligned(g_yz_xx_yz_x, g_yz_xx_yz_y, g_yz_xx_yz_z, g_z_0_y_0_y_xx_z_x, g_z_0_y_0_y_xx_z_y, g_z_0_y_0_y_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xx_z_x[i] = 4.0 * g_yz_xx_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_z_y[i] = 4.0 * g_yz_xx_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_z_z[i] = 4.0 * g_yz_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1197-1200)

    #pragma omp simd aligned(g_yz_xy_xy_x, g_yz_xy_xy_y, g_yz_xy_xy_z, g_z_0_y_0_y_xy_x_x, g_z_0_y_0_y_xy_x_y, g_z_0_y_0_y_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xy_x_x[i] = 4.0 * g_yz_xy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_x_y[i] = 4.0 * g_yz_xy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_x_z[i] = 4.0 * g_yz_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1200-1203)

    #pragma omp simd aligned(g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_yz_xy_yy_x, g_yz_xy_yy_y, g_yz_xy_yy_z, g_z_0_y_0_y_xy_y_x, g_z_0_y_0_y_xy_y_y, g_z_0_y_0_y_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xy_y_x[i] = -2.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_yz_xy_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_y_y[i] = -2.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_yz_xy_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_y_z[i] = -2.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_yz_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1203-1206)

    #pragma omp simd aligned(g_yz_xy_yz_x, g_yz_xy_yz_y, g_yz_xy_yz_z, g_z_0_y_0_y_xy_z_x, g_z_0_y_0_y_xy_z_y, g_z_0_y_0_y_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xy_z_x[i] = 4.0 * g_yz_xy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_z_y[i] = 4.0 * g_yz_xy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_z_z[i] = 4.0 * g_yz_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1206-1209)

    #pragma omp simd aligned(g_yz_xz_xy_x, g_yz_xz_xy_y, g_yz_xz_xy_z, g_z_0_y_0_y_xz_x_x, g_z_0_y_0_y_xz_x_y, g_z_0_y_0_y_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xz_x_x[i] = 4.0 * g_yz_xz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_x_y[i] = 4.0 * g_yz_xz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_x_z[i] = 4.0 * g_yz_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1209-1212)

    #pragma omp simd aligned(g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_yz_xz_yy_x, g_yz_xz_yy_y, g_yz_xz_yy_z, g_z_0_y_0_y_xz_y_x, g_z_0_y_0_y_xz_y_y, g_z_0_y_0_y_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xz_y_x[i] = -2.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_yz_xz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_y_y[i] = -2.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_yz_xz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_y_z[i] = -2.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_yz_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1212-1215)

    #pragma omp simd aligned(g_yz_xz_yz_x, g_yz_xz_yz_y, g_yz_xz_yz_z, g_z_0_y_0_y_xz_z_x, g_z_0_y_0_y_xz_z_y, g_z_0_y_0_y_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xz_z_x[i] = 4.0 * g_yz_xz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_z_y[i] = 4.0 * g_yz_xz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_z_z[i] = 4.0 * g_yz_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1215-1218)

    #pragma omp simd aligned(g_yz_yy_xy_x, g_yz_yy_xy_y, g_yz_yy_xy_z, g_z_0_y_0_y_yy_x_x, g_z_0_y_0_y_yy_x_y, g_z_0_y_0_y_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_yy_x_x[i] = 4.0 * g_yz_yy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_x_y[i] = 4.0 * g_yz_yy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_x_z[i] = 4.0 * g_yz_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1218-1221)

    #pragma omp simd aligned(g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z, g_yz_yy_yy_x, g_yz_yy_yy_y, g_yz_yy_yy_z, g_z_0_y_0_y_yy_y_x, g_z_0_y_0_y_yy_y_y, g_z_0_y_0_y_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_yy_y_x[i] = -2.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_yz_yy_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_y_y[i] = -2.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_yz_yy_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_y_z[i] = -2.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_yz_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1221-1224)

    #pragma omp simd aligned(g_yz_yy_yz_x, g_yz_yy_yz_y, g_yz_yy_yz_z, g_z_0_y_0_y_yy_z_x, g_z_0_y_0_y_yy_z_y, g_z_0_y_0_y_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_yy_z_x[i] = 4.0 * g_yz_yy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_z_y[i] = 4.0 * g_yz_yy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_z_z[i] = 4.0 * g_yz_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1224-1227)

    #pragma omp simd aligned(g_yz_yz_xy_x, g_yz_yz_xy_y, g_yz_yz_xy_z, g_z_0_y_0_y_yz_x_x, g_z_0_y_0_y_yz_x_y, g_z_0_y_0_y_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_yz_x_x[i] = 4.0 * g_yz_yz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_x_y[i] = 4.0 * g_yz_yz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_x_z[i] = 4.0 * g_yz_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1227-1230)

    #pragma omp simd aligned(g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_yz_yz_yy_x, g_yz_yz_yy_y, g_yz_yz_yy_z, g_z_0_y_0_y_yz_y_x, g_z_0_y_0_y_yz_y_y, g_z_0_y_0_y_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_yz_y_x[i] = -2.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_yz_yz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_y_y[i] = -2.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_yz_yz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_y_z[i] = -2.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_yz_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1230-1233)

    #pragma omp simd aligned(g_yz_yz_yz_x, g_yz_yz_yz_y, g_yz_yz_yz_z, g_z_0_y_0_y_yz_z_x, g_z_0_y_0_y_yz_z_y, g_z_0_y_0_y_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_yz_z_x[i] = 4.0 * g_yz_yz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_z_y[i] = 4.0 * g_yz_yz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_z_z[i] = 4.0 * g_yz_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1233-1236)

    #pragma omp simd aligned(g_yz_zz_xy_x, g_yz_zz_xy_y, g_yz_zz_xy_z, g_z_0_y_0_y_zz_x_x, g_z_0_y_0_y_zz_x_y, g_z_0_y_0_y_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_zz_x_x[i] = 4.0 * g_yz_zz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_x_y[i] = 4.0 * g_yz_zz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_x_z[i] = 4.0 * g_yz_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1236-1239)

    #pragma omp simd aligned(g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z, g_yz_zz_yy_x, g_yz_zz_yy_y, g_yz_zz_yy_z, g_z_0_y_0_y_zz_y_x, g_z_0_y_0_y_zz_y_y, g_z_0_y_0_y_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_zz_y_x[i] = -2.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_yz_zz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_y_y[i] = -2.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_yz_zz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_y_z[i] = -2.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_yz_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1239-1242)

    #pragma omp simd aligned(g_yz_zz_yz_x, g_yz_zz_yz_y, g_yz_zz_yz_z, g_z_0_y_0_y_zz_z_x, g_z_0_y_0_y_zz_z_y, g_z_0_y_0_y_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_zz_z_x[i] = 4.0 * g_yz_zz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_z_y[i] = 4.0 * g_yz_zz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_z_z[i] = 4.0 * g_yz_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1242-1245)

    #pragma omp simd aligned(g_0_xx_xy_x, g_0_xx_xy_y, g_0_xx_xy_z, g_z_0_y_0_z_xx_x_x, g_z_0_y_0_z_xx_x_y, g_z_0_y_0_z_xx_x_z, g_zz_xx_xy_x, g_zz_xx_xy_y, g_zz_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xx_x_x[i] = -2.0 * g_0_xx_xy_x[i] * c_exps[i] + 4.0 * g_zz_xx_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_x_y[i] = -2.0 * g_0_xx_xy_y[i] * c_exps[i] + 4.0 * g_zz_xx_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_x_z[i] = -2.0 * g_0_xx_xy_z[i] * c_exps[i] + 4.0 * g_zz_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1245-1248)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_0_xx_yy_x, g_0_xx_yy_y, g_0_xx_yy_z, g_z_0_y_0_z_xx_y_x, g_z_0_y_0_z_xx_y_y, g_z_0_y_0_z_xx_y_z, g_zz_xx_0_x, g_zz_xx_0_y, g_zz_xx_0_z, g_zz_xx_yy_x, g_zz_xx_yy_y, g_zz_xx_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xx_y_x[i] = g_0_xx_0_x[i] - 2.0 * g_0_xx_yy_x[i] * c_exps[i] - 2.0 * g_zz_xx_0_x[i] * a_exp + 4.0 * g_zz_xx_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_y_y[i] = g_0_xx_0_y[i] - 2.0 * g_0_xx_yy_y[i] * c_exps[i] - 2.0 * g_zz_xx_0_y[i] * a_exp + 4.0 * g_zz_xx_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_y_z[i] = g_0_xx_0_z[i] - 2.0 * g_0_xx_yy_z[i] * c_exps[i] - 2.0 * g_zz_xx_0_z[i] * a_exp + 4.0 * g_zz_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1248-1251)

    #pragma omp simd aligned(g_0_xx_yz_x, g_0_xx_yz_y, g_0_xx_yz_z, g_z_0_y_0_z_xx_z_x, g_z_0_y_0_z_xx_z_y, g_z_0_y_0_z_xx_z_z, g_zz_xx_yz_x, g_zz_xx_yz_y, g_zz_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xx_z_x[i] = -2.0 * g_0_xx_yz_x[i] * c_exps[i] + 4.0 * g_zz_xx_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_z_y[i] = -2.0 * g_0_xx_yz_y[i] * c_exps[i] + 4.0 * g_zz_xx_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_z_z[i] = -2.0 * g_0_xx_yz_z[i] * c_exps[i] + 4.0 * g_zz_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1251-1254)

    #pragma omp simd aligned(g_0_xy_xy_x, g_0_xy_xy_y, g_0_xy_xy_z, g_z_0_y_0_z_xy_x_x, g_z_0_y_0_z_xy_x_y, g_z_0_y_0_z_xy_x_z, g_zz_xy_xy_x, g_zz_xy_xy_y, g_zz_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xy_x_x[i] = -2.0 * g_0_xy_xy_x[i] * c_exps[i] + 4.0 * g_zz_xy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_x_y[i] = -2.0 * g_0_xy_xy_y[i] * c_exps[i] + 4.0 * g_zz_xy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_x_z[i] = -2.0 * g_0_xy_xy_z[i] * c_exps[i] + 4.0 * g_zz_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1254-1257)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_0_xy_yy_x, g_0_xy_yy_y, g_0_xy_yy_z, g_z_0_y_0_z_xy_y_x, g_z_0_y_0_z_xy_y_y, g_z_0_y_0_z_xy_y_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z, g_zz_xy_yy_x, g_zz_xy_yy_y, g_zz_xy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xy_y_x[i] = g_0_xy_0_x[i] - 2.0 * g_0_xy_yy_x[i] * c_exps[i] - 2.0 * g_zz_xy_0_x[i] * a_exp + 4.0 * g_zz_xy_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_y_y[i] = g_0_xy_0_y[i] - 2.0 * g_0_xy_yy_y[i] * c_exps[i] - 2.0 * g_zz_xy_0_y[i] * a_exp + 4.0 * g_zz_xy_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_y_z[i] = g_0_xy_0_z[i] - 2.0 * g_0_xy_yy_z[i] * c_exps[i] - 2.0 * g_zz_xy_0_z[i] * a_exp + 4.0 * g_zz_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1257-1260)

    #pragma omp simd aligned(g_0_xy_yz_x, g_0_xy_yz_y, g_0_xy_yz_z, g_z_0_y_0_z_xy_z_x, g_z_0_y_0_z_xy_z_y, g_z_0_y_0_z_xy_z_z, g_zz_xy_yz_x, g_zz_xy_yz_y, g_zz_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xy_z_x[i] = -2.0 * g_0_xy_yz_x[i] * c_exps[i] + 4.0 * g_zz_xy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_z_y[i] = -2.0 * g_0_xy_yz_y[i] * c_exps[i] + 4.0 * g_zz_xy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_z_z[i] = -2.0 * g_0_xy_yz_z[i] * c_exps[i] + 4.0 * g_zz_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1260-1263)

    #pragma omp simd aligned(g_0_xz_xy_x, g_0_xz_xy_y, g_0_xz_xy_z, g_z_0_y_0_z_xz_x_x, g_z_0_y_0_z_xz_x_y, g_z_0_y_0_z_xz_x_z, g_zz_xz_xy_x, g_zz_xz_xy_y, g_zz_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xz_x_x[i] = -2.0 * g_0_xz_xy_x[i] * c_exps[i] + 4.0 * g_zz_xz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_x_y[i] = -2.0 * g_0_xz_xy_y[i] * c_exps[i] + 4.0 * g_zz_xz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_x_z[i] = -2.0 * g_0_xz_xy_z[i] * c_exps[i] + 4.0 * g_zz_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1263-1266)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_0_xz_yy_x, g_0_xz_yy_y, g_0_xz_yy_z, g_z_0_y_0_z_xz_y_x, g_z_0_y_0_z_xz_y_y, g_z_0_y_0_z_xz_y_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z, g_zz_xz_yy_x, g_zz_xz_yy_y, g_zz_xz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xz_y_x[i] = g_0_xz_0_x[i] - 2.0 * g_0_xz_yy_x[i] * c_exps[i] - 2.0 * g_zz_xz_0_x[i] * a_exp + 4.0 * g_zz_xz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_y_y[i] = g_0_xz_0_y[i] - 2.0 * g_0_xz_yy_y[i] * c_exps[i] - 2.0 * g_zz_xz_0_y[i] * a_exp + 4.0 * g_zz_xz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_y_z[i] = g_0_xz_0_z[i] - 2.0 * g_0_xz_yy_z[i] * c_exps[i] - 2.0 * g_zz_xz_0_z[i] * a_exp + 4.0 * g_zz_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1266-1269)

    #pragma omp simd aligned(g_0_xz_yz_x, g_0_xz_yz_y, g_0_xz_yz_z, g_z_0_y_0_z_xz_z_x, g_z_0_y_0_z_xz_z_y, g_z_0_y_0_z_xz_z_z, g_zz_xz_yz_x, g_zz_xz_yz_y, g_zz_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xz_z_x[i] = -2.0 * g_0_xz_yz_x[i] * c_exps[i] + 4.0 * g_zz_xz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_z_y[i] = -2.0 * g_0_xz_yz_y[i] * c_exps[i] + 4.0 * g_zz_xz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_z_z[i] = -2.0 * g_0_xz_yz_z[i] * c_exps[i] + 4.0 * g_zz_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1269-1272)

    #pragma omp simd aligned(g_0_yy_xy_x, g_0_yy_xy_y, g_0_yy_xy_z, g_z_0_y_0_z_yy_x_x, g_z_0_y_0_z_yy_x_y, g_z_0_y_0_z_yy_x_z, g_zz_yy_xy_x, g_zz_yy_xy_y, g_zz_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_yy_x_x[i] = -2.0 * g_0_yy_xy_x[i] * c_exps[i] + 4.0 * g_zz_yy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_x_y[i] = -2.0 * g_0_yy_xy_y[i] * c_exps[i] + 4.0 * g_zz_yy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_x_z[i] = -2.0 * g_0_yy_xy_z[i] * c_exps[i] + 4.0 * g_zz_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1272-1275)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_0_yy_yy_x, g_0_yy_yy_y, g_0_yy_yy_z, g_z_0_y_0_z_yy_y_x, g_z_0_y_0_z_yy_y_y, g_z_0_y_0_z_yy_y_z, g_zz_yy_0_x, g_zz_yy_0_y, g_zz_yy_0_z, g_zz_yy_yy_x, g_zz_yy_yy_y, g_zz_yy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_yy_y_x[i] = g_0_yy_0_x[i] - 2.0 * g_0_yy_yy_x[i] * c_exps[i] - 2.0 * g_zz_yy_0_x[i] * a_exp + 4.0 * g_zz_yy_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_y_y[i] = g_0_yy_0_y[i] - 2.0 * g_0_yy_yy_y[i] * c_exps[i] - 2.0 * g_zz_yy_0_y[i] * a_exp + 4.0 * g_zz_yy_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_y_z[i] = g_0_yy_0_z[i] - 2.0 * g_0_yy_yy_z[i] * c_exps[i] - 2.0 * g_zz_yy_0_z[i] * a_exp + 4.0 * g_zz_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1275-1278)

    #pragma omp simd aligned(g_0_yy_yz_x, g_0_yy_yz_y, g_0_yy_yz_z, g_z_0_y_0_z_yy_z_x, g_z_0_y_0_z_yy_z_y, g_z_0_y_0_z_yy_z_z, g_zz_yy_yz_x, g_zz_yy_yz_y, g_zz_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_yy_z_x[i] = -2.0 * g_0_yy_yz_x[i] * c_exps[i] + 4.0 * g_zz_yy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_z_y[i] = -2.0 * g_0_yy_yz_y[i] * c_exps[i] + 4.0 * g_zz_yy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_z_z[i] = -2.0 * g_0_yy_yz_z[i] * c_exps[i] + 4.0 * g_zz_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1278-1281)

    #pragma omp simd aligned(g_0_yz_xy_x, g_0_yz_xy_y, g_0_yz_xy_z, g_z_0_y_0_z_yz_x_x, g_z_0_y_0_z_yz_x_y, g_z_0_y_0_z_yz_x_z, g_zz_yz_xy_x, g_zz_yz_xy_y, g_zz_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_yz_x_x[i] = -2.0 * g_0_yz_xy_x[i] * c_exps[i] + 4.0 * g_zz_yz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_x_y[i] = -2.0 * g_0_yz_xy_y[i] * c_exps[i] + 4.0 * g_zz_yz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_x_z[i] = -2.0 * g_0_yz_xy_z[i] * c_exps[i] + 4.0 * g_zz_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1281-1284)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_0_yz_yy_x, g_0_yz_yy_y, g_0_yz_yy_z, g_z_0_y_0_z_yz_y_x, g_z_0_y_0_z_yz_y_y, g_z_0_y_0_z_yz_y_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z, g_zz_yz_yy_x, g_zz_yz_yy_y, g_zz_yz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_yz_y_x[i] = g_0_yz_0_x[i] - 2.0 * g_0_yz_yy_x[i] * c_exps[i] - 2.0 * g_zz_yz_0_x[i] * a_exp + 4.0 * g_zz_yz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_y_y[i] = g_0_yz_0_y[i] - 2.0 * g_0_yz_yy_y[i] * c_exps[i] - 2.0 * g_zz_yz_0_y[i] * a_exp + 4.0 * g_zz_yz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_y_z[i] = g_0_yz_0_z[i] - 2.0 * g_0_yz_yy_z[i] * c_exps[i] - 2.0 * g_zz_yz_0_z[i] * a_exp + 4.0 * g_zz_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1284-1287)

    #pragma omp simd aligned(g_0_yz_yz_x, g_0_yz_yz_y, g_0_yz_yz_z, g_z_0_y_0_z_yz_z_x, g_z_0_y_0_z_yz_z_y, g_z_0_y_0_z_yz_z_z, g_zz_yz_yz_x, g_zz_yz_yz_y, g_zz_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_yz_z_x[i] = -2.0 * g_0_yz_yz_x[i] * c_exps[i] + 4.0 * g_zz_yz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_z_y[i] = -2.0 * g_0_yz_yz_y[i] * c_exps[i] + 4.0 * g_zz_yz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_z_z[i] = -2.0 * g_0_yz_yz_z[i] * c_exps[i] + 4.0 * g_zz_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1287-1290)

    #pragma omp simd aligned(g_0_zz_xy_x, g_0_zz_xy_y, g_0_zz_xy_z, g_z_0_y_0_z_zz_x_x, g_z_0_y_0_z_zz_x_y, g_z_0_y_0_z_zz_x_z, g_zz_zz_xy_x, g_zz_zz_xy_y, g_zz_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_zz_x_x[i] = -2.0 * g_0_zz_xy_x[i] * c_exps[i] + 4.0 * g_zz_zz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_x_y[i] = -2.0 * g_0_zz_xy_y[i] * c_exps[i] + 4.0 * g_zz_zz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_x_z[i] = -2.0 * g_0_zz_xy_z[i] * c_exps[i] + 4.0 * g_zz_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1290-1293)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_0_zz_yy_x, g_0_zz_yy_y, g_0_zz_yy_z, g_z_0_y_0_z_zz_y_x, g_z_0_y_0_z_zz_y_y, g_z_0_y_0_z_zz_y_z, g_zz_zz_0_x, g_zz_zz_0_y, g_zz_zz_0_z, g_zz_zz_yy_x, g_zz_zz_yy_y, g_zz_zz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_zz_y_x[i] = g_0_zz_0_x[i] - 2.0 * g_0_zz_yy_x[i] * c_exps[i] - 2.0 * g_zz_zz_0_x[i] * a_exp + 4.0 * g_zz_zz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_y_y[i] = g_0_zz_0_y[i] - 2.0 * g_0_zz_yy_y[i] * c_exps[i] - 2.0 * g_zz_zz_0_y[i] * a_exp + 4.0 * g_zz_zz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_y_z[i] = g_0_zz_0_z[i] - 2.0 * g_0_zz_yy_z[i] * c_exps[i] - 2.0 * g_zz_zz_0_z[i] * a_exp + 4.0 * g_zz_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1293-1296)

    #pragma omp simd aligned(g_0_zz_yz_x, g_0_zz_yz_y, g_0_zz_yz_z, g_z_0_y_0_z_zz_z_x, g_z_0_y_0_z_zz_z_y, g_z_0_y_0_z_zz_z_z, g_zz_zz_yz_x, g_zz_zz_yz_y, g_zz_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_zz_z_x[i] = -2.0 * g_0_zz_yz_x[i] * c_exps[i] + 4.0 * g_zz_zz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_z_y[i] = -2.0 * g_0_zz_yz_y[i] * c_exps[i] + 4.0 * g_zz_zz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_z_z[i] = -2.0 * g_0_zz_yz_z[i] * c_exps[i] + 4.0 * g_zz_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1296-1299)

    #pragma omp simd aligned(g_xz_xx_xz_x, g_xz_xx_xz_y, g_xz_xx_xz_z, g_z_0_z_0_x_xx_x_x, g_z_0_z_0_x_xx_x_y, g_z_0_z_0_x_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xx_x_x[i] = 4.0 * g_xz_xx_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_x_y[i] = 4.0 * g_xz_xx_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_x_z[i] = 4.0 * g_xz_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1299-1302)

    #pragma omp simd aligned(g_xz_xx_yz_x, g_xz_xx_yz_y, g_xz_xx_yz_z, g_z_0_z_0_x_xx_y_x, g_z_0_z_0_x_xx_y_y, g_z_0_z_0_x_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xx_y_x[i] = 4.0 * g_xz_xx_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_y_y[i] = 4.0 * g_xz_xx_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_y_z[i] = 4.0 * g_xz_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1302-1305)

    #pragma omp simd aligned(g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_xz_xx_zz_x, g_xz_xx_zz_y, g_xz_xx_zz_z, g_z_0_z_0_x_xx_z_x, g_z_0_z_0_x_xx_z_y, g_z_0_z_0_x_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xx_z_x[i] = -2.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xz_xx_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_z_y[i] = -2.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xz_xx_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_z_z[i] = -2.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xz_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1305-1308)

    #pragma omp simd aligned(g_xz_xy_xz_x, g_xz_xy_xz_y, g_xz_xy_xz_z, g_z_0_z_0_x_xy_x_x, g_z_0_z_0_x_xy_x_y, g_z_0_z_0_x_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xy_x_x[i] = 4.0 * g_xz_xy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_x_y[i] = 4.0 * g_xz_xy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_x_z[i] = 4.0 * g_xz_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1308-1311)

    #pragma omp simd aligned(g_xz_xy_yz_x, g_xz_xy_yz_y, g_xz_xy_yz_z, g_z_0_z_0_x_xy_y_x, g_z_0_z_0_x_xy_y_y, g_z_0_z_0_x_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xy_y_x[i] = 4.0 * g_xz_xy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_y_y[i] = 4.0 * g_xz_xy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_y_z[i] = 4.0 * g_xz_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1311-1314)

    #pragma omp simd aligned(g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_xz_xy_zz_x, g_xz_xy_zz_y, g_xz_xy_zz_z, g_z_0_z_0_x_xy_z_x, g_z_0_z_0_x_xy_z_y, g_z_0_z_0_x_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xy_z_x[i] = -2.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xz_xy_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_z_y[i] = -2.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xz_xy_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_z_z[i] = -2.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xz_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1314-1317)

    #pragma omp simd aligned(g_xz_xz_xz_x, g_xz_xz_xz_y, g_xz_xz_xz_z, g_z_0_z_0_x_xz_x_x, g_z_0_z_0_x_xz_x_y, g_z_0_z_0_x_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xz_x_x[i] = 4.0 * g_xz_xz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_x_y[i] = 4.0 * g_xz_xz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_x_z[i] = 4.0 * g_xz_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1317-1320)

    #pragma omp simd aligned(g_xz_xz_yz_x, g_xz_xz_yz_y, g_xz_xz_yz_z, g_z_0_z_0_x_xz_y_x, g_z_0_z_0_x_xz_y_y, g_z_0_z_0_x_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xz_y_x[i] = 4.0 * g_xz_xz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_y_y[i] = 4.0 * g_xz_xz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_y_z[i] = 4.0 * g_xz_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1320-1323)

    #pragma omp simd aligned(g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_xz_xz_zz_x, g_xz_xz_zz_y, g_xz_xz_zz_z, g_z_0_z_0_x_xz_z_x, g_z_0_z_0_x_xz_z_y, g_z_0_z_0_x_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xz_z_x[i] = -2.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xz_xz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_z_y[i] = -2.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xz_xz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_z_z[i] = -2.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xz_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1323-1326)

    #pragma omp simd aligned(g_xz_yy_xz_x, g_xz_yy_xz_y, g_xz_yy_xz_z, g_z_0_z_0_x_yy_x_x, g_z_0_z_0_x_yy_x_y, g_z_0_z_0_x_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_yy_x_x[i] = 4.0 * g_xz_yy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_x_y[i] = 4.0 * g_xz_yy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_x_z[i] = 4.0 * g_xz_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1326-1329)

    #pragma omp simd aligned(g_xz_yy_yz_x, g_xz_yy_yz_y, g_xz_yy_yz_z, g_z_0_z_0_x_yy_y_x, g_z_0_z_0_x_yy_y_y, g_z_0_z_0_x_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_yy_y_x[i] = 4.0 * g_xz_yy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_y_y[i] = 4.0 * g_xz_yy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_y_z[i] = 4.0 * g_xz_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1329-1332)

    #pragma omp simd aligned(g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_xz_yy_zz_x, g_xz_yy_zz_y, g_xz_yy_zz_z, g_z_0_z_0_x_yy_z_x, g_z_0_z_0_x_yy_z_y, g_z_0_z_0_x_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_yy_z_x[i] = -2.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xz_yy_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_z_y[i] = -2.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xz_yy_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_z_z[i] = -2.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xz_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1332-1335)

    #pragma omp simd aligned(g_xz_yz_xz_x, g_xz_yz_xz_y, g_xz_yz_xz_z, g_z_0_z_0_x_yz_x_x, g_z_0_z_0_x_yz_x_y, g_z_0_z_0_x_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_yz_x_x[i] = 4.0 * g_xz_yz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_x_y[i] = 4.0 * g_xz_yz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_x_z[i] = 4.0 * g_xz_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1335-1338)

    #pragma omp simd aligned(g_xz_yz_yz_x, g_xz_yz_yz_y, g_xz_yz_yz_z, g_z_0_z_0_x_yz_y_x, g_z_0_z_0_x_yz_y_y, g_z_0_z_0_x_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_yz_y_x[i] = 4.0 * g_xz_yz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_y_y[i] = 4.0 * g_xz_yz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_y_z[i] = 4.0 * g_xz_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1338-1341)

    #pragma omp simd aligned(g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_xz_yz_zz_x, g_xz_yz_zz_y, g_xz_yz_zz_z, g_z_0_z_0_x_yz_z_x, g_z_0_z_0_x_yz_z_y, g_z_0_z_0_x_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_yz_z_x[i] = -2.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xz_yz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_z_y[i] = -2.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xz_yz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_z_z[i] = -2.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xz_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1341-1344)

    #pragma omp simd aligned(g_xz_zz_xz_x, g_xz_zz_xz_y, g_xz_zz_xz_z, g_z_0_z_0_x_zz_x_x, g_z_0_z_0_x_zz_x_y, g_z_0_z_0_x_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_zz_x_x[i] = 4.0 * g_xz_zz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_x_y[i] = 4.0 * g_xz_zz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_x_z[i] = 4.0 * g_xz_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1344-1347)

    #pragma omp simd aligned(g_xz_zz_yz_x, g_xz_zz_yz_y, g_xz_zz_yz_z, g_z_0_z_0_x_zz_y_x, g_z_0_z_0_x_zz_y_y, g_z_0_z_0_x_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_zz_y_x[i] = 4.0 * g_xz_zz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_y_y[i] = 4.0 * g_xz_zz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_y_z[i] = 4.0 * g_xz_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1347-1350)

    #pragma omp simd aligned(g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_xz_zz_zz_x, g_xz_zz_zz_y, g_xz_zz_zz_z, g_z_0_z_0_x_zz_z_x, g_z_0_z_0_x_zz_z_y, g_z_0_z_0_x_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_zz_z_x[i] = -2.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xz_zz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_z_y[i] = -2.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xz_zz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_z_z[i] = -2.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xz_zz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1350-1353)

    #pragma omp simd aligned(g_yz_xx_xz_x, g_yz_xx_xz_y, g_yz_xx_xz_z, g_z_0_z_0_y_xx_x_x, g_z_0_z_0_y_xx_x_y, g_z_0_z_0_y_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xx_x_x[i] = 4.0 * g_yz_xx_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_x_y[i] = 4.0 * g_yz_xx_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_x_z[i] = 4.0 * g_yz_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1353-1356)

    #pragma omp simd aligned(g_yz_xx_yz_x, g_yz_xx_yz_y, g_yz_xx_yz_z, g_z_0_z_0_y_xx_y_x, g_z_0_z_0_y_xx_y_y, g_z_0_z_0_y_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xx_y_x[i] = 4.0 * g_yz_xx_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_y_y[i] = 4.0 * g_yz_xx_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_y_z[i] = 4.0 * g_yz_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1356-1359)

    #pragma omp simd aligned(g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z, g_yz_xx_zz_x, g_yz_xx_zz_y, g_yz_xx_zz_z, g_z_0_z_0_y_xx_z_x, g_z_0_z_0_y_xx_z_y, g_z_0_z_0_y_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xx_z_x[i] = -2.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_yz_xx_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_z_y[i] = -2.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_yz_xx_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_z_z[i] = -2.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_yz_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1359-1362)

    #pragma omp simd aligned(g_yz_xy_xz_x, g_yz_xy_xz_y, g_yz_xy_xz_z, g_z_0_z_0_y_xy_x_x, g_z_0_z_0_y_xy_x_y, g_z_0_z_0_y_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xy_x_x[i] = 4.0 * g_yz_xy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_x_y[i] = 4.0 * g_yz_xy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_x_z[i] = 4.0 * g_yz_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1362-1365)

    #pragma omp simd aligned(g_yz_xy_yz_x, g_yz_xy_yz_y, g_yz_xy_yz_z, g_z_0_z_0_y_xy_y_x, g_z_0_z_0_y_xy_y_y, g_z_0_z_0_y_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xy_y_x[i] = 4.0 * g_yz_xy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_y_y[i] = 4.0 * g_yz_xy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_y_z[i] = 4.0 * g_yz_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1365-1368)

    #pragma omp simd aligned(g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_yz_xy_zz_x, g_yz_xy_zz_y, g_yz_xy_zz_z, g_z_0_z_0_y_xy_z_x, g_z_0_z_0_y_xy_z_y, g_z_0_z_0_y_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xy_z_x[i] = -2.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_yz_xy_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_z_y[i] = -2.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_yz_xy_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_z_z[i] = -2.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_yz_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1368-1371)

    #pragma omp simd aligned(g_yz_xz_xz_x, g_yz_xz_xz_y, g_yz_xz_xz_z, g_z_0_z_0_y_xz_x_x, g_z_0_z_0_y_xz_x_y, g_z_0_z_0_y_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xz_x_x[i] = 4.0 * g_yz_xz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_x_y[i] = 4.0 * g_yz_xz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_x_z[i] = 4.0 * g_yz_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1371-1374)

    #pragma omp simd aligned(g_yz_xz_yz_x, g_yz_xz_yz_y, g_yz_xz_yz_z, g_z_0_z_0_y_xz_y_x, g_z_0_z_0_y_xz_y_y, g_z_0_z_0_y_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xz_y_x[i] = 4.0 * g_yz_xz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_y_y[i] = 4.0 * g_yz_xz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_y_z[i] = 4.0 * g_yz_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1374-1377)

    #pragma omp simd aligned(g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_yz_xz_zz_x, g_yz_xz_zz_y, g_yz_xz_zz_z, g_z_0_z_0_y_xz_z_x, g_z_0_z_0_y_xz_z_y, g_z_0_z_0_y_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xz_z_x[i] = -2.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_yz_xz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_z_y[i] = -2.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_yz_xz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_z_z[i] = -2.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_yz_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1377-1380)

    #pragma omp simd aligned(g_yz_yy_xz_x, g_yz_yy_xz_y, g_yz_yy_xz_z, g_z_0_z_0_y_yy_x_x, g_z_0_z_0_y_yy_x_y, g_z_0_z_0_y_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_yy_x_x[i] = 4.0 * g_yz_yy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_x_y[i] = 4.0 * g_yz_yy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_x_z[i] = 4.0 * g_yz_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1380-1383)

    #pragma omp simd aligned(g_yz_yy_yz_x, g_yz_yy_yz_y, g_yz_yy_yz_z, g_z_0_z_0_y_yy_y_x, g_z_0_z_0_y_yy_y_y, g_z_0_z_0_y_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_yy_y_x[i] = 4.0 * g_yz_yy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_y_y[i] = 4.0 * g_yz_yy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_y_z[i] = 4.0 * g_yz_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1383-1386)

    #pragma omp simd aligned(g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z, g_yz_yy_zz_x, g_yz_yy_zz_y, g_yz_yy_zz_z, g_z_0_z_0_y_yy_z_x, g_z_0_z_0_y_yy_z_y, g_z_0_z_0_y_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_yy_z_x[i] = -2.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_yz_yy_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_z_y[i] = -2.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_yz_yy_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_z_z[i] = -2.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_yz_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1386-1389)

    #pragma omp simd aligned(g_yz_yz_xz_x, g_yz_yz_xz_y, g_yz_yz_xz_z, g_z_0_z_0_y_yz_x_x, g_z_0_z_0_y_yz_x_y, g_z_0_z_0_y_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_yz_x_x[i] = 4.0 * g_yz_yz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_x_y[i] = 4.0 * g_yz_yz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_x_z[i] = 4.0 * g_yz_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1389-1392)

    #pragma omp simd aligned(g_yz_yz_yz_x, g_yz_yz_yz_y, g_yz_yz_yz_z, g_z_0_z_0_y_yz_y_x, g_z_0_z_0_y_yz_y_y, g_z_0_z_0_y_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_yz_y_x[i] = 4.0 * g_yz_yz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_y_y[i] = 4.0 * g_yz_yz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_y_z[i] = 4.0 * g_yz_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1392-1395)

    #pragma omp simd aligned(g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_yz_yz_zz_x, g_yz_yz_zz_y, g_yz_yz_zz_z, g_z_0_z_0_y_yz_z_x, g_z_0_z_0_y_yz_z_y, g_z_0_z_0_y_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_yz_z_x[i] = -2.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_yz_yz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_z_y[i] = -2.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_yz_yz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_z_z[i] = -2.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_yz_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1395-1398)

    #pragma omp simd aligned(g_yz_zz_xz_x, g_yz_zz_xz_y, g_yz_zz_xz_z, g_z_0_z_0_y_zz_x_x, g_z_0_z_0_y_zz_x_y, g_z_0_z_0_y_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_zz_x_x[i] = 4.0 * g_yz_zz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_x_y[i] = 4.0 * g_yz_zz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_x_z[i] = 4.0 * g_yz_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1398-1401)

    #pragma omp simd aligned(g_yz_zz_yz_x, g_yz_zz_yz_y, g_yz_zz_yz_z, g_z_0_z_0_y_zz_y_x, g_z_0_z_0_y_zz_y_y, g_z_0_z_0_y_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_zz_y_x[i] = 4.0 * g_yz_zz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_y_y[i] = 4.0 * g_yz_zz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_y_z[i] = 4.0 * g_yz_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1401-1404)

    #pragma omp simd aligned(g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z, g_yz_zz_zz_x, g_yz_zz_zz_y, g_yz_zz_zz_z, g_z_0_z_0_y_zz_z_x, g_z_0_z_0_y_zz_z_y, g_z_0_z_0_y_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_zz_z_x[i] = -2.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_yz_zz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_z_y[i] = -2.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_yz_zz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_z_z[i] = -2.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_yz_zz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1404-1407)

    #pragma omp simd aligned(g_0_xx_xz_x, g_0_xx_xz_y, g_0_xx_xz_z, g_z_0_z_0_z_xx_x_x, g_z_0_z_0_z_xx_x_y, g_z_0_z_0_z_xx_x_z, g_zz_xx_xz_x, g_zz_xx_xz_y, g_zz_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xx_x_x[i] = -2.0 * g_0_xx_xz_x[i] * c_exps[i] + 4.0 * g_zz_xx_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_x_y[i] = -2.0 * g_0_xx_xz_y[i] * c_exps[i] + 4.0 * g_zz_xx_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_x_z[i] = -2.0 * g_0_xx_xz_z[i] * c_exps[i] + 4.0 * g_zz_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1407-1410)

    #pragma omp simd aligned(g_0_xx_yz_x, g_0_xx_yz_y, g_0_xx_yz_z, g_z_0_z_0_z_xx_y_x, g_z_0_z_0_z_xx_y_y, g_z_0_z_0_z_xx_y_z, g_zz_xx_yz_x, g_zz_xx_yz_y, g_zz_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xx_y_x[i] = -2.0 * g_0_xx_yz_x[i] * c_exps[i] + 4.0 * g_zz_xx_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_y_y[i] = -2.0 * g_0_xx_yz_y[i] * c_exps[i] + 4.0 * g_zz_xx_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_y_z[i] = -2.0 * g_0_xx_yz_z[i] * c_exps[i] + 4.0 * g_zz_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1410-1413)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_0_xx_zz_x, g_0_xx_zz_y, g_0_xx_zz_z, g_z_0_z_0_z_xx_z_x, g_z_0_z_0_z_xx_z_y, g_z_0_z_0_z_xx_z_z, g_zz_xx_0_x, g_zz_xx_0_y, g_zz_xx_0_z, g_zz_xx_zz_x, g_zz_xx_zz_y, g_zz_xx_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xx_z_x[i] = g_0_xx_0_x[i] - 2.0 * g_0_xx_zz_x[i] * c_exps[i] - 2.0 * g_zz_xx_0_x[i] * a_exp + 4.0 * g_zz_xx_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_z_y[i] = g_0_xx_0_y[i] - 2.0 * g_0_xx_zz_y[i] * c_exps[i] - 2.0 * g_zz_xx_0_y[i] * a_exp + 4.0 * g_zz_xx_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_z_z[i] = g_0_xx_0_z[i] - 2.0 * g_0_xx_zz_z[i] * c_exps[i] - 2.0 * g_zz_xx_0_z[i] * a_exp + 4.0 * g_zz_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1413-1416)

    #pragma omp simd aligned(g_0_xy_xz_x, g_0_xy_xz_y, g_0_xy_xz_z, g_z_0_z_0_z_xy_x_x, g_z_0_z_0_z_xy_x_y, g_z_0_z_0_z_xy_x_z, g_zz_xy_xz_x, g_zz_xy_xz_y, g_zz_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xy_x_x[i] = -2.0 * g_0_xy_xz_x[i] * c_exps[i] + 4.0 * g_zz_xy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_x_y[i] = -2.0 * g_0_xy_xz_y[i] * c_exps[i] + 4.0 * g_zz_xy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_x_z[i] = -2.0 * g_0_xy_xz_z[i] * c_exps[i] + 4.0 * g_zz_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1416-1419)

    #pragma omp simd aligned(g_0_xy_yz_x, g_0_xy_yz_y, g_0_xy_yz_z, g_z_0_z_0_z_xy_y_x, g_z_0_z_0_z_xy_y_y, g_z_0_z_0_z_xy_y_z, g_zz_xy_yz_x, g_zz_xy_yz_y, g_zz_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xy_y_x[i] = -2.0 * g_0_xy_yz_x[i] * c_exps[i] + 4.0 * g_zz_xy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_y_y[i] = -2.0 * g_0_xy_yz_y[i] * c_exps[i] + 4.0 * g_zz_xy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_y_z[i] = -2.0 * g_0_xy_yz_z[i] * c_exps[i] + 4.0 * g_zz_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1419-1422)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_0_xy_zz_x, g_0_xy_zz_y, g_0_xy_zz_z, g_z_0_z_0_z_xy_z_x, g_z_0_z_0_z_xy_z_y, g_z_0_z_0_z_xy_z_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z, g_zz_xy_zz_x, g_zz_xy_zz_y, g_zz_xy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xy_z_x[i] = g_0_xy_0_x[i] - 2.0 * g_0_xy_zz_x[i] * c_exps[i] - 2.0 * g_zz_xy_0_x[i] * a_exp + 4.0 * g_zz_xy_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_z_y[i] = g_0_xy_0_y[i] - 2.0 * g_0_xy_zz_y[i] * c_exps[i] - 2.0 * g_zz_xy_0_y[i] * a_exp + 4.0 * g_zz_xy_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_z_z[i] = g_0_xy_0_z[i] - 2.0 * g_0_xy_zz_z[i] * c_exps[i] - 2.0 * g_zz_xy_0_z[i] * a_exp + 4.0 * g_zz_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1422-1425)

    #pragma omp simd aligned(g_0_xz_xz_x, g_0_xz_xz_y, g_0_xz_xz_z, g_z_0_z_0_z_xz_x_x, g_z_0_z_0_z_xz_x_y, g_z_0_z_0_z_xz_x_z, g_zz_xz_xz_x, g_zz_xz_xz_y, g_zz_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xz_x_x[i] = -2.0 * g_0_xz_xz_x[i] * c_exps[i] + 4.0 * g_zz_xz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_x_y[i] = -2.0 * g_0_xz_xz_y[i] * c_exps[i] + 4.0 * g_zz_xz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_x_z[i] = -2.0 * g_0_xz_xz_z[i] * c_exps[i] + 4.0 * g_zz_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1425-1428)

    #pragma omp simd aligned(g_0_xz_yz_x, g_0_xz_yz_y, g_0_xz_yz_z, g_z_0_z_0_z_xz_y_x, g_z_0_z_0_z_xz_y_y, g_z_0_z_0_z_xz_y_z, g_zz_xz_yz_x, g_zz_xz_yz_y, g_zz_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xz_y_x[i] = -2.0 * g_0_xz_yz_x[i] * c_exps[i] + 4.0 * g_zz_xz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_y_y[i] = -2.0 * g_0_xz_yz_y[i] * c_exps[i] + 4.0 * g_zz_xz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_y_z[i] = -2.0 * g_0_xz_yz_z[i] * c_exps[i] + 4.0 * g_zz_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1428-1431)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_0_xz_zz_x, g_0_xz_zz_y, g_0_xz_zz_z, g_z_0_z_0_z_xz_z_x, g_z_0_z_0_z_xz_z_y, g_z_0_z_0_z_xz_z_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z, g_zz_xz_zz_x, g_zz_xz_zz_y, g_zz_xz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xz_z_x[i] = g_0_xz_0_x[i] - 2.0 * g_0_xz_zz_x[i] * c_exps[i] - 2.0 * g_zz_xz_0_x[i] * a_exp + 4.0 * g_zz_xz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_z_y[i] = g_0_xz_0_y[i] - 2.0 * g_0_xz_zz_y[i] * c_exps[i] - 2.0 * g_zz_xz_0_y[i] * a_exp + 4.0 * g_zz_xz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_z_z[i] = g_0_xz_0_z[i] - 2.0 * g_0_xz_zz_z[i] * c_exps[i] - 2.0 * g_zz_xz_0_z[i] * a_exp + 4.0 * g_zz_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1431-1434)

    #pragma omp simd aligned(g_0_yy_xz_x, g_0_yy_xz_y, g_0_yy_xz_z, g_z_0_z_0_z_yy_x_x, g_z_0_z_0_z_yy_x_y, g_z_0_z_0_z_yy_x_z, g_zz_yy_xz_x, g_zz_yy_xz_y, g_zz_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_yy_x_x[i] = -2.0 * g_0_yy_xz_x[i] * c_exps[i] + 4.0 * g_zz_yy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_x_y[i] = -2.0 * g_0_yy_xz_y[i] * c_exps[i] + 4.0 * g_zz_yy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_x_z[i] = -2.0 * g_0_yy_xz_z[i] * c_exps[i] + 4.0 * g_zz_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1434-1437)

    #pragma omp simd aligned(g_0_yy_yz_x, g_0_yy_yz_y, g_0_yy_yz_z, g_z_0_z_0_z_yy_y_x, g_z_0_z_0_z_yy_y_y, g_z_0_z_0_z_yy_y_z, g_zz_yy_yz_x, g_zz_yy_yz_y, g_zz_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_yy_y_x[i] = -2.0 * g_0_yy_yz_x[i] * c_exps[i] + 4.0 * g_zz_yy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_y_y[i] = -2.0 * g_0_yy_yz_y[i] * c_exps[i] + 4.0 * g_zz_yy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_y_z[i] = -2.0 * g_0_yy_yz_z[i] * c_exps[i] + 4.0 * g_zz_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1437-1440)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_0_yy_zz_x, g_0_yy_zz_y, g_0_yy_zz_z, g_z_0_z_0_z_yy_z_x, g_z_0_z_0_z_yy_z_y, g_z_0_z_0_z_yy_z_z, g_zz_yy_0_x, g_zz_yy_0_y, g_zz_yy_0_z, g_zz_yy_zz_x, g_zz_yy_zz_y, g_zz_yy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_yy_z_x[i] = g_0_yy_0_x[i] - 2.0 * g_0_yy_zz_x[i] * c_exps[i] - 2.0 * g_zz_yy_0_x[i] * a_exp + 4.0 * g_zz_yy_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_z_y[i] = g_0_yy_0_y[i] - 2.0 * g_0_yy_zz_y[i] * c_exps[i] - 2.0 * g_zz_yy_0_y[i] * a_exp + 4.0 * g_zz_yy_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_z_z[i] = g_0_yy_0_z[i] - 2.0 * g_0_yy_zz_z[i] * c_exps[i] - 2.0 * g_zz_yy_0_z[i] * a_exp + 4.0 * g_zz_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1440-1443)

    #pragma omp simd aligned(g_0_yz_xz_x, g_0_yz_xz_y, g_0_yz_xz_z, g_z_0_z_0_z_yz_x_x, g_z_0_z_0_z_yz_x_y, g_z_0_z_0_z_yz_x_z, g_zz_yz_xz_x, g_zz_yz_xz_y, g_zz_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_yz_x_x[i] = -2.0 * g_0_yz_xz_x[i] * c_exps[i] + 4.0 * g_zz_yz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_x_y[i] = -2.0 * g_0_yz_xz_y[i] * c_exps[i] + 4.0 * g_zz_yz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_x_z[i] = -2.0 * g_0_yz_xz_z[i] * c_exps[i] + 4.0 * g_zz_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1443-1446)

    #pragma omp simd aligned(g_0_yz_yz_x, g_0_yz_yz_y, g_0_yz_yz_z, g_z_0_z_0_z_yz_y_x, g_z_0_z_0_z_yz_y_y, g_z_0_z_0_z_yz_y_z, g_zz_yz_yz_x, g_zz_yz_yz_y, g_zz_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_yz_y_x[i] = -2.0 * g_0_yz_yz_x[i] * c_exps[i] + 4.0 * g_zz_yz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_y_y[i] = -2.0 * g_0_yz_yz_y[i] * c_exps[i] + 4.0 * g_zz_yz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_y_z[i] = -2.0 * g_0_yz_yz_z[i] * c_exps[i] + 4.0 * g_zz_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1446-1449)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_0_yz_zz_x, g_0_yz_zz_y, g_0_yz_zz_z, g_z_0_z_0_z_yz_z_x, g_z_0_z_0_z_yz_z_y, g_z_0_z_0_z_yz_z_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z, g_zz_yz_zz_x, g_zz_yz_zz_y, g_zz_yz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_yz_z_x[i] = g_0_yz_0_x[i] - 2.0 * g_0_yz_zz_x[i] * c_exps[i] - 2.0 * g_zz_yz_0_x[i] * a_exp + 4.0 * g_zz_yz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_z_y[i] = g_0_yz_0_y[i] - 2.0 * g_0_yz_zz_y[i] * c_exps[i] - 2.0 * g_zz_yz_0_y[i] * a_exp + 4.0 * g_zz_yz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_z_z[i] = g_0_yz_0_z[i] - 2.0 * g_0_yz_zz_z[i] * c_exps[i] - 2.0 * g_zz_yz_0_z[i] * a_exp + 4.0 * g_zz_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1449-1452)

    #pragma omp simd aligned(g_0_zz_xz_x, g_0_zz_xz_y, g_0_zz_xz_z, g_z_0_z_0_z_zz_x_x, g_z_0_z_0_z_zz_x_y, g_z_0_z_0_z_zz_x_z, g_zz_zz_xz_x, g_zz_zz_xz_y, g_zz_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_zz_x_x[i] = -2.0 * g_0_zz_xz_x[i] * c_exps[i] + 4.0 * g_zz_zz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_x_y[i] = -2.0 * g_0_zz_xz_y[i] * c_exps[i] + 4.0 * g_zz_zz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_x_z[i] = -2.0 * g_0_zz_xz_z[i] * c_exps[i] + 4.0 * g_zz_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1452-1455)

    #pragma omp simd aligned(g_0_zz_yz_x, g_0_zz_yz_y, g_0_zz_yz_z, g_z_0_z_0_z_zz_y_x, g_z_0_z_0_z_zz_y_y, g_z_0_z_0_z_zz_y_z, g_zz_zz_yz_x, g_zz_zz_yz_y, g_zz_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_zz_y_x[i] = -2.0 * g_0_zz_yz_x[i] * c_exps[i] + 4.0 * g_zz_zz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_y_y[i] = -2.0 * g_0_zz_yz_y[i] * c_exps[i] + 4.0 * g_zz_zz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_y_z[i] = -2.0 * g_0_zz_yz_z[i] * c_exps[i] + 4.0 * g_zz_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (1455-1458)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_0_zz_zz_x, g_0_zz_zz_y, g_0_zz_zz_z, g_z_0_z_0_z_zz_z_x, g_z_0_z_0_z_zz_z_y, g_z_0_z_0_z_zz_z_z, g_zz_zz_0_x, g_zz_zz_0_y, g_zz_zz_0_z, g_zz_zz_zz_x, g_zz_zz_zz_y, g_zz_zz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_zz_z_x[i] = g_0_zz_0_x[i] - 2.0 * g_0_zz_zz_x[i] * c_exps[i] - 2.0 * g_zz_zz_0_x[i] * a_exp + 4.0 * g_zz_zz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_z_y[i] = g_0_zz_0_y[i] - 2.0 * g_0_zz_zz_y[i] * c_exps[i] - 2.0 * g_zz_zz_0_y[i] * a_exp + 4.0 * g_zz_zz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_z_z[i] = g_0_zz_0_z[i] - 2.0 * g_0_zz_zz_z[i] * c_exps[i] - 2.0 * g_zz_zz_0_z[i] * a_exp + 4.0 * g_zz_zz_zz_z[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

