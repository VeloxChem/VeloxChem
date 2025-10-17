#include "ElectronRepulsionGeom0010ContrRecXXHF.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxhf(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxhf,
                                            const size_t idx_xxgf,
                                            const size_t idx_geom_10_xxgf,
                                            const size_t idx_geom_10_xxgg,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSGF

            const auto gf_off = idx_xxgf + (i * bcomps + j) * 150;

            auto g_xxxx_xxx = cbuffer.data(gf_off + 0);

            auto g_xxxx_xxy = cbuffer.data(gf_off + 1);

            auto g_xxxx_xxz = cbuffer.data(gf_off + 2);

            auto g_xxxx_xyy = cbuffer.data(gf_off + 3);

            auto g_xxxx_xyz = cbuffer.data(gf_off + 4);

            auto g_xxxx_xzz = cbuffer.data(gf_off + 5);

            auto g_xxxx_yyy = cbuffer.data(gf_off + 6);

            auto g_xxxx_yyz = cbuffer.data(gf_off + 7);

            auto g_xxxx_yzz = cbuffer.data(gf_off + 8);

            auto g_xxxx_zzz = cbuffer.data(gf_off + 9);

            auto g_yyyy_xxx = cbuffer.data(gf_off + 100);

            auto g_yyyy_xxy = cbuffer.data(gf_off + 101);

            auto g_yyyy_xxz = cbuffer.data(gf_off + 102);

            auto g_yyyy_xyy = cbuffer.data(gf_off + 103);

            auto g_yyyy_xyz = cbuffer.data(gf_off + 104);

            auto g_yyyy_xzz = cbuffer.data(gf_off + 105);

            auto g_yyyy_yyy = cbuffer.data(gf_off + 106);

            auto g_yyyy_yyz = cbuffer.data(gf_off + 107);

            auto g_yyyy_yzz = cbuffer.data(gf_off + 108);

            auto g_yyyy_zzz = cbuffer.data(gf_off + 109);

            auto g_zzzz_xxx = cbuffer.data(gf_off + 140);

            auto g_zzzz_xxy = cbuffer.data(gf_off + 141);

            auto g_zzzz_xxz = cbuffer.data(gf_off + 142);

            auto g_zzzz_xyy = cbuffer.data(gf_off + 143);

            auto g_zzzz_xyz = cbuffer.data(gf_off + 144);

            auto g_zzzz_xzz = cbuffer.data(gf_off + 145);

            auto g_zzzz_yyy = cbuffer.data(gf_off + 146);

            auto g_zzzz_yyz = cbuffer.data(gf_off + 147);

            auto g_zzzz_yzz = cbuffer.data(gf_off + 148);

            auto g_zzzz_zzz = cbuffer.data(gf_off + 149);

            /// Set up components of auxilary buffer : SSGF

            const auto gf_geom_10_off = idx_geom_10_xxgf + (i * bcomps + j) * 150;

            auto g_x_0_xxxx_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxx_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxx_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxx_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxx_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxy_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxy_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxy_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxy_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxy_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxy_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxy_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxy_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxy_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxy_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxyy_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxyy_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxyy_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxyy_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxyy_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxyy_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxyy_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxyy_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxyy_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxyy_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxyz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxyz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxyz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxyz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxyz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxyz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxyz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxyz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxyz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxyz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xyyy_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xyyy_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xyyy_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xyyy_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xyyy_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xyyy_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xyyy_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xyyy_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xyyy_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xyyy_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xyyz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xyyz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xyyz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xyyz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xyyz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xyyz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xyyz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xyyz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xyyz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xyyz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xyzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xyzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xyzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xyzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xyzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xyzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xyzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xyzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xyzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xyzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xzzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xzzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xzzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xzzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xzzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xzzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xzzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xzzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xzzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xzzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_yyyy_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_yyyy_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_yyyy_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_yyyy_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_yyyy_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_yyyy_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_yyyy_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_yyyy_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_yyyy_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_yyyy_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_yyyz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_yyyz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_yyyz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_yyyz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_yyyz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_yyyz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_yyyz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_yyyz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_yyyz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_yyyz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_yyzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_yyzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_yyzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_yyzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_yyzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_yyzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_yyzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_yyzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_yyzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_yyzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_yzzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_yzzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_yzzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_yzzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_yzzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_yzzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_yzzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_yzzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_yzzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_yzzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_zzzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_zzzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_zzzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_zzzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_zzzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_zzzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_zzzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_zzzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_zzzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_zzzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_y_0_xxxx_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 2);

            auto g_y_0_xxxx_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 3);

            auto g_y_0_xxxx_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 4);

            auto g_y_0_xxxx_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 5);

            auto g_y_0_xxxx_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 6);

            auto g_y_0_xxxx_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 7);

            auto g_y_0_xxxx_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 8);

            auto g_y_0_xxxx_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 9);

            auto g_y_0_xxxy_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 10);

            auto g_y_0_xxxy_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 11);

            auto g_y_0_xxxy_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 12);

            auto g_y_0_xxxy_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 13);

            auto g_y_0_xxxy_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 14);

            auto g_y_0_xxxy_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 15);

            auto g_y_0_xxxy_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 16);

            auto g_y_0_xxxy_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 17);

            auto g_y_0_xxxy_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 18);

            auto g_y_0_xxxy_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 19);

            auto g_y_0_xxxz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 20);

            auto g_y_0_xxxz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 21);

            auto g_y_0_xxxz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 22);

            auto g_y_0_xxxz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 23);

            auto g_y_0_xxxz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 24);

            auto g_y_0_xxxz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 25);

            auto g_y_0_xxxz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 26);

            auto g_y_0_xxxz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 27);

            auto g_y_0_xxxz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 28);

            auto g_y_0_xxxz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 29);

            auto g_y_0_xxyy_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 30);

            auto g_y_0_xxyy_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 31);

            auto g_y_0_xxyy_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 32);

            auto g_y_0_xxyy_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 33);

            auto g_y_0_xxyy_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 34);

            auto g_y_0_xxyy_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 35);

            auto g_y_0_xxyy_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 36);

            auto g_y_0_xxyy_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 37);

            auto g_y_0_xxyy_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 38);

            auto g_y_0_xxyy_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 39);

            auto g_y_0_xxyz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 40);

            auto g_y_0_xxyz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 41);

            auto g_y_0_xxyz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 42);

            auto g_y_0_xxyz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 43);

            auto g_y_0_xxyz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 44);

            auto g_y_0_xxyz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 45);

            auto g_y_0_xxyz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 46);

            auto g_y_0_xxyz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 47);

            auto g_y_0_xxyz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 48);

            auto g_y_0_xxyz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 49);

            auto g_y_0_xxzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 50);

            auto g_y_0_xxzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 51);

            auto g_y_0_xxzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 52);

            auto g_y_0_xxzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 53);

            auto g_y_0_xxzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 54);

            auto g_y_0_xxzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 55);

            auto g_y_0_xxzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 56);

            auto g_y_0_xxzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 57);

            auto g_y_0_xxzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 58);

            auto g_y_0_xxzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 59);

            auto g_y_0_xyyy_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 60);

            auto g_y_0_xyyy_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 61);

            auto g_y_0_xyyy_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 62);

            auto g_y_0_xyyy_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 63);

            auto g_y_0_xyyy_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 64);

            auto g_y_0_xyyy_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 65);

            auto g_y_0_xyyy_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 66);

            auto g_y_0_xyyy_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 67);

            auto g_y_0_xyyy_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 68);

            auto g_y_0_xyyy_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 69);

            auto g_y_0_xyyz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 70);

            auto g_y_0_xyyz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 71);

            auto g_y_0_xyyz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 72);

            auto g_y_0_xyyz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 73);

            auto g_y_0_xyyz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 74);

            auto g_y_0_xyyz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 75);

            auto g_y_0_xyyz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 76);

            auto g_y_0_xyyz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 77);

            auto g_y_0_xyyz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 78);

            auto g_y_0_xyyz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 79);

            auto g_y_0_xyzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 80);

            auto g_y_0_xyzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 81);

            auto g_y_0_xyzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 82);

            auto g_y_0_xyzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 83);

            auto g_y_0_xyzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 84);

            auto g_y_0_xyzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 85);

            auto g_y_0_xyzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 86);

            auto g_y_0_xyzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 87);

            auto g_y_0_xyzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 88);

            auto g_y_0_xyzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 89);

            auto g_y_0_xzzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 90);

            auto g_y_0_xzzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 91);

            auto g_y_0_xzzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 92);

            auto g_y_0_xzzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 93);

            auto g_y_0_xzzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 94);

            auto g_y_0_xzzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 95);

            auto g_y_0_xzzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 96);

            auto g_y_0_xzzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 97);

            auto g_y_0_xzzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 98);

            auto g_y_0_xzzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 99);

            auto g_y_0_yyyy_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 100);

            auto g_y_0_yyyy_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 101);

            auto g_y_0_yyyy_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 102);

            auto g_y_0_yyyy_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 103);

            auto g_y_0_yyyy_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 104);

            auto g_y_0_yyyy_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 105);

            auto g_y_0_yyyy_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 106);

            auto g_y_0_yyyy_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 107);

            auto g_y_0_yyyy_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 108);

            auto g_y_0_yyyy_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 109);

            auto g_y_0_yyyz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 110);

            auto g_y_0_yyyz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 111);

            auto g_y_0_yyyz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 112);

            auto g_y_0_yyyz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 113);

            auto g_y_0_yyyz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 114);

            auto g_y_0_yyyz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 115);

            auto g_y_0_yyyz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 116);

            auto g_y_0_yyyz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 117);

            auto g_y_0_yyyz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 118);

            auto g_y_0_yyyz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 119);

            auto g_y_0_yyzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 120);

            auto g_y_0_yyzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 121);

            auto g_y_0_yyzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 122);

            auto g_y_0_yyzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 123);

            auto g_y_0_yyzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 124);

            auto g_y_0_yyzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 125);

            auto g_y_0_yyzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 126);

            auto g_y_0_yyzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 127);

            auto g_y_0_yyzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 128);

            auto g_y_0_yyzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 129);

            auto g_y_0_yzzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 130);

            auto g_y_0_yzzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 131);

            auto g_y_0_yzzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 132);

            auto g_y_0_yzzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 133);

            auto g_y_0_yzzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 134);

            auto g_y_0_yzzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 135);

            auto g_y_0_yzzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 136);

            auto g_y_0_yzzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 137);

            auto g_y_0_yzzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 138);

            auto g_y_0_yzzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 139);

            auto g_y_0_zzzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 140);

            auto g_y_0_zzzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 141);

            auto g_y_0_zzzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 142);

            auto g_y_0_zzzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 143);

            auto g_y_0_zzzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 144);

            auto g_y_0_zzzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 145);

            auto g_y_0_zzzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 146);

            auto g_y_0_zzzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 147);

            auto g_y_0_zzzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 148);

            auto g_y_0_zzzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 149);

            auto g_z_0_xxxx_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 2);

            auto g_z_0_xxxx_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 3);

            auto g_z_0_xxxx_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 4);

            auto g_z_0_xxxx_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 5);

            auto g_z_0_xxxx_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 6);

            auto g_z_0_xxxx_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 7);

            auto g_z_0_xxxx_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 8);

            auto g_z_0_xxxx_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 9);

            auto g_z_0_xxxy_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 10);

            auto g_z_0_xxxy_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 11);

            auto g_z_0_xxxy_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 12);

            auto g_z_0_xxxy_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 13);

            auto g_z_0_xxxy_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 14);

            auto g_z_0_xxxy_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 15);

            auto g_z_0_xxxy_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 16);

            auto g_z_0_xxxy_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 17);

            auto g_z_0_xxxy_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 18);

            auto g_z_0_xxxy_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 19);

            auto g_z_0_xxxz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 20);

            auto g_z_0_xxxz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 21);

            auto g_z_0_xxxz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 22);

            auto g_z_0_xxxz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 23);

            auto g_z_0_xxxz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 24);

            auto g_z_0_xxxz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 25);

            auto g_z_0_xxxz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 26);

            auto g_z_0_xxxz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 27);

            auto g_z_0_xxxz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 28);

            auto g_z_0_xxxz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 29);

            auto g_z_0_xxyy_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 30);

            auto g_z_0_xxyy_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 31);

            auto g_z_0_xxyy_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 32);

            auto g_z_0_xxyy_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 33);

            auto g_z_0_xxyy_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 34);

            auto g_z_0_xxyy_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 35);

            auto g_z_0_xxyy_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 36);

            auto g_z_0_xxyy_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 37);

            auto g_z_0_xxyy_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 38);

            auto g_z_0_xxyy_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 39);

            auto g_z_0_xxyz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 40);

            auto g_z_0_xxyz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 41);

            auto g_z_0_xxyz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 42);

            auto g_z_0_xxyz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 43);

            auto g_z_0_xxyz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 44);

            auto g_z_0_xxyz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 45);

            auto g_z_0_xxyz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 46);

            auto g_z_0_xxyz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 47);

            auto g_z_0_xxyz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 48);

            auto g_z_0_xxyz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 49);

            auto g_z_0_xxzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 50);

            auto g_z_0_xxzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 51);

            auto g_z_0_xxzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 52);

            auto g_z_0_xxzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 53);

            auto g_z_0_xxzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 54);

            auto g_z_0_xxzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 55);

            auto g_z_0_xxzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 56);

            auto g_z_0_xxzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 57);

            auto g_z_0_xxzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 58);

            auto g_z_0_xxzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 59);

            auto g_z_0_xyyy_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 60);

            auto g_z_0_xyyy_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 61);

            auto g_z_0_xyyy_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 62);

            auto g_z_0_xyyy_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 63);

            auto g_z_0_xyyy_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 64);

            auto g_z_0_xyyy_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 65);

            auto g_z_0_xyyy_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 66);

            auto g_z_0_xyyy_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 67);

            auto g_z_0_xyyy_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 68);

            auto g_z_0_xyyy_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 69);

            auto g_z_0_xyyz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 70);

            auto g_z_0_xyyz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 71);

            auto g_z_0_xyyz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 72);

            auto g_z_0_xyyz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 73);

            auto g_z_0_xyyz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 74);

            auto g_z_0_xyyz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 75);

            auto g_z_0_xyyz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 76);

            auto g_z_0_xyyz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 77);

            auto g_z_0_xyyz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 78);

            auto g_z_0_xyyz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 79);

            auto g_z_0_xyzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 80);

            auto g_z_0_xyzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 81);

            auto g_z_0_xyzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 82);

            auto g_z_0_xyzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 83);

            auto g_z_0_xyzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 84);

            auto g_z_0_xyzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 85);

            auto g_z_0_xyzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 86);

            auto g_z_0_xyzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 87);

            auto g_z_0_xyzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 88);

            auto g_z_0_xyzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 89);

            auto g_z_0_xzzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 90);

            auto g_z_0_xzzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 91);

            auto g_z_0_xzzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 92);

            auto g_z_0_xzzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 93);

            auto g_z_0_xzzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 94);

            auto g_z_0_xzzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 95);

            auto g_z_0_xzzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 96);

            auto g_z_0_xzzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 97);

            auto g_z_0_xzzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 98);

            auto g_z_0_xzzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 99);

            auto g_z_0_yyyy_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 100);

            auto g_z_0_yyyy_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 101);

            auto g_z_0_yyyy_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 102);

            auto g_z_0_yyyy_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 103);

            auto g_z_0_yyyy_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 104);

            auto g_z_0_yyyy_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 105);

            auto g_z_0_yyyy_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 106);

            auto g_z_0_yyyy_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 107);

            auto g_z_0_yyyy_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 108);

            auto g_z_0_yyyy_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 109);

            auto g_z_0_yyyz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 110);

            auto g_z_0_yyyz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 111);

            auto g_z_0_yyyz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 112);

            auto g_z_0_yyyz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 113);

            auto g_z_0_yyyz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 114);

            auto g_z_0_yyyz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 115);

            auto g_z_0_yyyz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 116);

            auto g_z_0_yyyz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 117);

            auto g_z_0_yyyz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 118);

            auto g_z_0_yyyz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 119);

            auto g_z_0_yyzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 120);

            auto g_z_0_yyzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 121);

            auto g_z_0_yyzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 122);

            auto g_z_0_yyzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 123);

            auto g_z_0_yyzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 124);

            auto g_z_0_yyzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 125);

            auto g_z_0_yyzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 126);

            auto g_z_0_yyzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 127);

            auto g_z_0_yyzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 128);

            auto g_z_0_yyzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 129);

            auto g_z_0_yzzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 130);

            auto g_z_0_yzzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 131);

            auto g_z_0_yzzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 132);

            auto g_z_0_yzzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 133);

            auto g_z_0_yzzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 134);

            auto g_z_0_yzzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 135);

            auto g_z_0_yzzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 136);

            auto g_z_0_yzzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 137);

            auto g_z_0_yzzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 138);

            auto g_z_0_yzzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 139);

            auto g_z_0_zzzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 140);

            auto g_z_0_zzzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 141);

            auto g_z_0_zzzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 142);

            auto g_z_0_zzzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 143);

            auto g_z_0_zzzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 144);

            auto g_z_0_zzzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 145);

            auto g_z_0_zzzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 146);

            auto g_z_0_zzzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 147);

            auto g_z_0_zzzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 148);

            auto g_z_0_zzzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 149);

            /// Set up components of auxilary buffer : SSGG

            const auto gg_geom_10_off = idx_geom_10_xxgg + (i * bcomps + j) * 225;

            auto g_x_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_y_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 2);

            auto g_y_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 3);

            auto g_y_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 4);

            auto g_y_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 5);

            auto g_y_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 6);

            auto g_y_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 7);

            auto g_y_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 8);

            auto g_y_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 9);

            auto g_y_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 15);

            auto g_y_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 16);

            auto g_y_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 17);

            auto g_y_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 18);

            auto g_y_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 19);

            auto g_y_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 20);

            auto g_y_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 21);

            auto g_y_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 22);

            auto g_y_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 23);

            auto g_y_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 24);

            auto g_y_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 30);

            auto g_y_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 31);

            auto g_y_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 32);

            auto g_y_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 33);

            auto g_y_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 34);

            auto g_y_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 35);

            auto g_y_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 36);

            auto g_y_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 37);

            auto g_y_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 38);

            auto g_y_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 39);

            auto g_y_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 45);

            auto g_y_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 46);

            auto g_y_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 47);

            auto g_y_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 48);

            auto g_y_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 49);

            auto g_y_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 50);

            auto g_y_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 51);

            auto g_y_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 52);

            auto g_y_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 53);

            auto g_y_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 54);

            auto g_y_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 60);

            auto g_y_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 61);

            auto g_y_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 62);

            auto g_y_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 63);

            auto g_y_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 64);

            auto g_y_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 65);

            auto g_y_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 66);

            auto g_y_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 67);

            auto g_y_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 68);

            auto g_y_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 69);

            auto g_y_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 75);

            auto g_y_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 76);

            auto g_y_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 77);

            auto g_y_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 78);

            auto g_y_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 79);

            auto g_y_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 80);

            auto g_y_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 81);

            auto g_y_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 82);

            auto g_y_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 83);

            auto g_y_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 84);

            auto g_y_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 90);

            auto g_y_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 91);

            auto g_y_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 92);

            auto g_y_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 93);

            auto g_y_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 94);

            auto g_y_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 95);

            auto g_y_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 96);

            auto g_y_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 97);

            auto g_y_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 98);

            auto g_y_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 99);

            auto g_y_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 105);

            auto g_y_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 106);

            auto g_y_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 107);

            auto g_y_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 108);

            auto g_y_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 109);

            auto g_y_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 110);

            auto g_y_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 111);

            auto g_y_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 112);

            auto g_y_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 113);

            auto g_y_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 114);

            auto g_y_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 120);

            auto g_y_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 121);

            auto g_y_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 122);

            auto g_y_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 123);

            auto g_y_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 124);

            auto g_y_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 125);

            auto g_y_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 126);

            auto g_y_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 127);

            auto g_y_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 128);

            auto g_y_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 129);

            auto g_y_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 135);

            auto g_y_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 136);

            auto g_y_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 137);

            auto g_y_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 138);

            auto g_y_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 139);

            auto g_y_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 140);

            auto g_y_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 141);

            auto g_y_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 142);

            auto g_y_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 143);

            auto g_y_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 144);

            auto g_y_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 150);

            auto g_y_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 151);

            auto g_y_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 152);

            auto g_y_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 153);

            auto g_y_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 154);

            auto g_y_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 155);

            auto g_y_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 156);

            auto g_y_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 157);

            auto g_y_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 158);

            auto g_y_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 159);

            auto g_y_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 160);

            auto g_y_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 161);

            auto g_y_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 162);

            auto g_y_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 163);

            auto g_y_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 164);

            auto g_y_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 165);

            auto g_y_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 166);

            auto g_y_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 167);

            auto g_y_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 168);

            auto g_y_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 169);

            auto g_y_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 170);

            auto g_y_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 171);

            auto g_y_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 172);

            auto g_y_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 173);

            auto g_y_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 174);

            auto g_y_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 176);

            auto g_y_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 177);

            auto g_y_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 178);

            auto g_y_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 179);

            auto g_y_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 180);

            auto g_y_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 181);

            auto g_y_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 182);

            auto g_y_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 183);

            auto g_y_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 184);

            auto g_y_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 185);

            auto g_y_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 186);

            auto g_y_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 187);

            auto g_y_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 188);

            auto g_y_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 189);

            auto g_y_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 191);

            auto g_y_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 192);

            auto g_y_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 193);

            auto g_y_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 194);

            auto g_y_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 195);

            auto g_y_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 196);

            auto g_y_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 197);

            auto g_y_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 198);

            auto g_y_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 199);

            auto g_y_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 200);

            auto g_y_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 201);

            auto g_y_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 202);

            auto g_y_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 203);

            auto g_y_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 204);

            auto g_y_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 206);

            auto g_y_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 207);

            auto g_y_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 208);

            auto g_y_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 209);

            auto g_y_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 210);

            auto g_y_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 211);

            auto g_y_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 212);

            auto g_y_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 213);

            auto g_y_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 214);

            auto g_y_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 215);

            auto g_y_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 216);

            auto g_y_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 217);

            auto g_y_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 218);

            auto g_y_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 219);

            auto g_y_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 221);

            auto g_y_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 222);

            auto g_y_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 223);

            auto g_y_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 224);

            auto g_z_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 2);

            auto g_z_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 3);

            auto g_z_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 4);

            auto g_z_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 5);

            auto g_z_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 6);

            auto g_z_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 7);

            auto g_z_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 8);

            auto g_z_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 9);

            auto g_z_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 15);

            auto g_z_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 16);

            auto g_z_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 17);

            auto g_z_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 18);

            auto g_z_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 19);

            auto g_z_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 20);

            auto g_z_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 21);

            auto g_z_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 22);

            auto g_z_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 23);

            auto g_z_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 24);

            auto g_z_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 30);

            auto g_z_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 31);

            auto g_z_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 32);

            auto g_z_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 33);

            auto g_z_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 34);

            auto g_z_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 35);

            auto g_z_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 36);

            auto g_z_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 37);

            auto g_z_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 38);

            auto g_z_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 39);

            auto g_z_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 45);

            auto g_z_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 46);

            auto g_z_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 47);

            auto g_z_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 48);

            auto g_z_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 49);

            auto g_z_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 50);

            auto g_z_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 51);

            auto g_z_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 52);

            auto g_z_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 53);

            auto g_z_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 54);

            auto g_z_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 60);

            auto g_z_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 61);

            auto g_z_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 62);

            auto g_z_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 63);

            auto g_z_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 64);

            auto g_z_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 65);

            auto g_z_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 66);

            auto g_z_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 67);

            auto g_z_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 68);

            auto g_z_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 69);

            auto g_z_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 75);

            auto g_z_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 76);

            auto g_z_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 77);

            auto g_z_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 78);

            auto g_z_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 79);

            auto g_z_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 80);

            auto g_z_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 81);

            auto g_z_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 82);

            auto g_z_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 83);

            auto g_z_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 84);

            auto g_z_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 90);

            auto g_z_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 91);

            auto g_z_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 92);

            auto g_z_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 93);

            auto g_z_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 94);

            auto g_z_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 95);

            auto g_z_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 96);

            auto g_z_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 97);

            auto g_z_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 98);

            auto g_z_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 99);

            auto g_z_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 105);

            auto g_z_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 106);

            auto g_z_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 107);

            auto g_z_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 108);

            auto g_z_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 109);

            auto g_z_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 110);

            auto g_z_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 111);

            auto g_z_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 112);

            auto g_z_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 113);

            auto g_z_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 114);

            auto g_z_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 120);

            auto g_z_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 121);

            auto g_z_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 122);

            auto g_z_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 123);

            auto g_z_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 124);

            auto g_z_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 125);

            auto g_z_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 126);

            auto g_z_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 127);

            auto g_z_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 128);

            auto g_z_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 129);

            auto g_z_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 135);

            auto g_z_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 136);

            auto g_z_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 137);

            auto g_z_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 138);

            auto g_z_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 139);

            auto g_z_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 140);

            auto g_z_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 141);

            auto g_z_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 142);

            auto g_z_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 143);

            auto g_z_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 144);

            auto g_z_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 150);

            auto g_z_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 151);

            auto g_z_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 152);

            auto g_z_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 153);

            auto g_z_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 154);

            auto g_z_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 155);

            auto g_z_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 156);

            auto g_z_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 157);

            auto g_z_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 158);

            auto g_z_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 159);

            auto g_z_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 160);

            auto g_z_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 161);

            auto g_z_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 162);

            auto g_z_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 163);

            auto g_z_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 165);

            auto g_z_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 166);

            auto g_z_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 167);

            auto g_z_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 168);

            auto g_z_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 169);

            auto g_z_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 170);

            auto g_z_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 171);

            auto g_z_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 172);

            auto g_z_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 173);

            auto g_z_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 174);

            auto g_z_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 175);

            auto g_z_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 176);

            auto g_z_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 177);

            auto g_z_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 178);

            auto g_z_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 180);

            auto g_z_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 181);

            auto g_z_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 182);

            auto g_z_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 183);

            auto g_z_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 184);

            auto g_z_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 185);

            auto g_z_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 186);

            auto g_z_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 187);

            auto g_z_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 188);

            auto g_z_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 189);

            auto g_z_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 190);

            auto g_z_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 191);

            auto g_z_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 192);

            auto g_z_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 193);

            auto g_z_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 195);

            auto g_z_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 196);

            auto g_z_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 197);

            auto g_z_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 198);

            auto g_z_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 199);

            auto g_z_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 200);

            auto g_z_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 201);

            auto g_z_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 202);

            auto g_z_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 203);

            auto g_z_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 204);

            auto g_z_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 205);

            auto g_z_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 206);

            auto g_z_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 207);

            auto g_z_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 208);

            auto g_z_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 210);

            auto g_z_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 211);

            auto g_z_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 212);

            auto g_z_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 213);

            auto g_z_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 214);

            auto g_z_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 215);

            auto g_z_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 216);

            auto g_z_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 217);

            auto g_z_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 218);

            auto g_z_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 219);

            auto g_z_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 220);

            auto g_z_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 221);

            auto g_z_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 222);

            auto g_z_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 223);

            auto g_z_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 224);

            /// set up bra offset for contr_buffer_xxhf

            const auto hf_geom_10_off = idx_geom_10_xxhf + (i * bcomps + j) * 210;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_x_0_xxxx_xxx, g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxy, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xyy, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_yyy, g_x_0_xxxx_yyz, g_x_0_xxxx_yzz, g_x_0_xxxx_zzz, g_x_0_xxxxx_xxx, g_x_0_xxxxx_xxy, g_x_0_xxxxx_xxz, g_x_0_xxxxx_xyy, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xzz, g_x_0_xxxxx_yyy, g_x_0_xxxxx_yyz, g_x_0_xxxxx_yzz, g_x_0_xxxxx_zzz, g_xxxx_xxx, g_xxxx_xxy, g_xxxx_xxz, g_xxxx_xyy, g_xxxx_xyz, g_xxxx_xzz, g_xxxx_yyy, g_xxxx_yyz, g_xxxx_yzz, g_xxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_xxx[k] = -g_xxxx_xxx[k] - g_x_0_xxxx_xxx[k] * cd_x[k] + g_x_0_xxxx_xxxx[k];

                g_x_0_xxxxx_xxy[k] = -g_xxxx_xxy[k] - g_x_0_xxxx_xxy[k] * cd_x[k] + g_x_0_xxxx_xxxy[k];

                g_x_0_xxxxx_xxz[k] = -g_xxxx_xxz[k] - g_x_0_xxxx_xxz[k] * cd_x[k] + g_x_0_xxxx_xxxz[k];

                g_x_0_xxxxx_xyy[k] = -g_xxxx_xyy[k] - g_x_0_xxxx_xyy[k] * cd_x[k] + g_x_0_xxxx_xxyy[k];

                g_x_0_xxxxx_xyz[k] = -g_xxxx_xyz[k] - g_x_0_xxxx_xyz[k] * cd_x[k] + g_x_0_xxxx_xxyz[k];

                g_x_0_xxxxx_xzz[k] = -g_xxxx_xzz[k] - g_x_0_xxxx_xzz[k] * cd_x[k] + g_x_0_xxxx_xxzz[k];

                g_x_0_xxxxx_yyy[k] = -g_xxxx_yyy[k] - g_x_0_xxxx_yyy[k] * cd_x[k] + g_x_0_xxxx_xyyy[k];

                g_x_0_xxxxx_yyz[k] = -g_xxxx_yyz[k] - g_x_0_xxxx_yyz[k] * cd_x[k] + g_x_0_xxxx_xyyz[k];

                g_x_0_xxxxx_yzz[k] = -g_xxxx_yzz[k] - g_x_0_xxxx_yzz[k] * cd_x[k] + g_x_0_xxxx_xyzz[k];

                g_x_0_xxxxx_zzz[k] = -g_xxxx_zzz[k] - g_x_0_xxxx_zzz[k] * cd_x[k] + g_x_0_xxxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_y, g_x_0_xxxx_xxx, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxy, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxz, g_x_0_xxxx_xyy, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xzz, g_x_0_xxxx_yyy, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_zzz, g_x_0_xxxxy_xxx, g_x_0_xxxxy_xxy, g_x_0_xxxxy_xxz, g_x_0_xxxxy_xyy, g_x_0_xxxxy_xyz, g_x_0_xxxxy_xzz, g_x_0_xxxxy_yyy, g_x_0_xxxxy_yyz, g_x_0_xxxxy_yzz, g_x_0_xxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_xxx[k] = -g_x_0_xxxx_xxx[k] * cd_y[k] + g_x_0_xxxx_xxxy[k];

                g_x_0_xxxxy_xxy[k] = -g_x_0_xxxx_xxy[k] * cd_y[k] + g_x_0_xxxx_xxyy[k];

                g_x_0_xxxxy_xxz[k] = -g_x_0_xxxx_xxz[k] * cd_y[k] + g_x_0_xxxx_xxyz[k];

                g_x_0_xxxxy_xyy[k] = -g_x_0_xxxx_xyy[k] * cd_y[k] + g_x_0_xxxx_xyyy[k];

                g_x_0_xxxxy_xyz[k] = -g_x_0_xxxx_xyz[k] * cd_y[k] + g_x_0_xxxx_xyyz[k];

                g_x_0_xxxxy_xzz[k] = -g_x_0_xxxx_xzz[k] * cd_y[k] + g_x_0_xxxx_xyzz[k];

                g_x_0_xxxxy_yyy[k] = -g_x_0_xxxx_yyy[k] * cd_y[k] + g_x_0_xxxx_yyyy[k];

                g_x_0_xxxxy_yyz[k] = -g_x_0_xxxx_yyz[k] * cd_y[k] + g_x_0_xxxx_yyyz[k];

                g_x_0_xxxxy_yzz[k] = -g_x_0_xxxx_yzz[k] * cd_y[k] + g_x_0_xxxx_yyzz[k];

                g_x_0_xxxxy_zzz[k] = -g_x_0_xxxx_zzz[k] * cd_y[k] + g_x_0_xxxx_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_z, g_x_0_xxxx_xxx, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxy, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xyy, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_yyy, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_zzz, g_x_0_xxxx_zzzz, g_x_0_xxxxz_xxx, g_x_0_xxxxz_xxy, g_x_0_xxxxz_xxz, g_x_0_xxxxz_xyy, g_x_0_xxxxz_xyz, g_x_0_xxxxz_xzz, g_x_0_xxxxz_yyy, g_x_0_xxxxz_yyz, g_x_0_xxxxz_yzz, g_x_0_xxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_xxx[k] = -g_x_0_xxxx_xxx[k] * cd_z[k] + g_x_0_xxxx_xxxz[k];

                g_x_0_xxxxz_xxy[k] = -g_x_0_xxxx_xxy[k] * cd_z[k] + g_x_0_xxxx_xxyz[k];

                g_x_0_xxxxz_xxz[k] = -g_x_0_xxxx_xxz[k] * cd_z[k] + g_x_0_xxxx_xxzz[k];

                g_x_0_xxxxz_xyy[k] = -g_x_0_xxxx_xyy[k] * cd_z[k] + g_x_0_xxxx_xyyz[k];

                g_x_0_xxxxz_xyz[k] = -g_x_0_xxxx_xyz[k] * cd_z[k] + g_x_0_xxxx_xyzz[k];

                g_x_0_xxxxz_xzz[k] = -g_x_0_xxxx_xzz[k] * cd_z[k] + g_x_0_xxxx_xzzz[k];

                g_x_0_xxxxz_yyy[k] = -g_x_0_xxxx_yyy[k] * cd_z[k] + g_x_0_xxxx_yyyz[k];

                g_x_0_xxxxz_yyz[k] = -g_x_0_xxxx_yyz[k] * cd_z[k] + g_x_0_xxxx_yyzz[k];

                g_x_0_xxxxz_yzz[k] = -g_x_0_xxxx_yzz[k] * cd_z[k] + g_x_0_xxxx_yzzz[k];

                g_x_0_xxxxz_zzz[k] = -g_x_0_xxxx_zzz[k] * cd_z[k] + g_x_0_xxxx_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 39);

            #pragma omp simd aligned(cd_y, g_x_0_xxxy_xxx, g_x_0_xxxy_xxxy, g_x_0_xxxy_xxy, g_x_0_xxxy_xxyy, g_x_0_xxxy_xxyz, g_x_0_xxxy_xxz, g_x_0_xxxy_xyy, g_x_0_xxxy_xyyy, g_x_0_xxxy_xyyz, g_x_0_xxxy_xyz, g_x_0_xxxy_xyzz, g_x_0_xxxy_xzz, g_x_0_xxxy_yyy, g_x_0_xxxy_yyyy, g_x_0_xxxy_yyyz, g_x_0_xxxy_yyz, g_x_0_xxxy_yyzz, g_x_0_xxxy_yzz, g_x_0_xxxy_yzzz, g_x_0_xxxy_zzz, g_x_0_xxxyy_xxx, g_x_0_xxxyy_xxy, g_x_0_xxxyy_xxz, g_x_0_xxxyy_xyy, g_x_0_xxxyy_xyz, g_x_0_xxxyy_xzz, g_x_0_xxxyy_yyy, g_x_0_xxxyy_yyz, g_x_0_xxxyy_yzz, g_x_0_xxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_xxx[k] = -g_x_0_xxxy_xxx[k] * cd_y[k] + g_x_0_xxxy_xxxy[k];

                g_x_0_xxxyy_xxy[k] = -g_x_0_xxxy_xxy[k] * cd_y[k] + g_x_0_xxxy_xxyy[k];

                g_x_0_xxxyy_xxz[k] = -g_x_0_xxxy_xxz[k] * cd_y[k] + g_x_0_xxxy_xxyz[k];

                g_x_0_xxxyy_xyy[k] = -g_x_0_xxxy_xyy[k] * cd_y[k] + g_x_0_xxxy_xyyy[k];

                g_x_0_xxxyy_xyz[k] = -g_x_0_xxxy_xyz[k] * cd_y[k] + g_x_0_xxxy_xyyz[k];

                g_x_0_xxxyy_xzz[k] = -g_x_0_xxxy_xzz[k] * cd_y[k] + g_x_0_xxxy_xyzz[k];

                g_x_0_xxxyy_yyy[k] = -g_x_0_xxxy_yyy[k] * cd_y[k] + g_x_0_xxxy_yyyy[k];

                g_x_0_xxxyy_yyz[k] = -g_x_0_xxxy_yyz[k] * cd_y[k] + g_x_0_xxxy_yyyz[k];

                g_x_0_xxxyy_yzz[k] = -g_x_0_xxxy_yzz[k] * cd_y[k] + g_x_0_xxxy_yyzz[k];

                g_x_0_xxxyy_zzz[k] = -g_x_0_xxxy_zzz[k] * cd_y[k] + g_x_0_xxxy_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 49);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyz_xxx, g_x_0_xxxyz_xxy, g_x_0_xxxyz_xxz, g_x_0_xxxyz_xyy, g_x_0_xxxyz_xyz, g_x_0_xxxyz_xzz, g_x_0_xxxyz_yyy, g_x_0_xxxyz_yyz, g_x_0_xxxyz_yzz, g_x_0_xxxyz_zzz, g_x_0_xxxz_xxx, g_x_0_xxxz_xxxy, g_x_0_xxxz_xxy, g_x_0_xxxz_xxyy, g_x_0_xxxz_xxyz, g_x_0_xxxz_xxz, g_x_0_xxxz_xyy, g_x_0_xxxz_xyyy, g_x_0_xxxz_xyyz, g_x_0_xxxz_xyz, g_x_0_xxxz_xyzz, g_x_0_xxxz_xzz, g_x_0_xxxz_yyy, g_x_0_xxxz_yyyy, g_x_0_xxxz_yyyz, g_x_0_xxxz_yyz, g_x_0_xxxz_yyzz, g_x_0_xxxz_yzz, g_x_0_xxxz_yzzz, g_x_0_xxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_xxx[k] = -g_x_0_xxxz_xxx[k] * cd_y[k] + g_x_0_xxxz_xxxy[k];

                g_x_0_xxxyz_xxy[k] = -g_x_0_xxxz_xxy[k] * cd_y[k] + g_x_0_xxxz_xxyy[k];

                g_x_0_xxxyz_xxz[k] = -g_x_0_xxxz_xxz[k] * cd_y[k] + g_x_0_xxxz_xxyz[k];

                g_x_0_xxxyz_xyy[k] = -g_x_0_xxxz_xyy[k] * cd_y[k] + g_x_0_xxxz_xyyy[k];

                g_x_0_xxxyz_xyz[k] = -g_x_0_xxxz_xyz[k] * cd_y[k] + g_x_0_xxxz_xyyz[k];

                g_x_0_xxxyz_xzz[k] = -g_x_0_xxxz_xzz[k] * cd_y[k] + g_x_0_xxxz_xyzz[k];

                g_x_0_xxxyz_yyy[k] = -g_x_0_xxxz_yyy[k] * cd_y[k] + g_x_0_xxxz_yyyy[k];

                g_x_0_xxxyz_yyz[k] = -g_x_0_xxxz_yyz[k] * cd_y[k] + g_x_0_xxxz_yyyz[k];

                g_x_0_xxxyz_yzz[k] = -g_x_0_xxxz_yzz[k] * cd_y[k] + g_x_0_xxxz_yyzz[k];

                g_x_0_xxxyz_zzz[k] = -g_x_0_xxxz_zzz[k] * cd_y[k] + g_x_0_xxxz_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_z, g_x_0_xxxz_xxx, g_x_0_xxxz_xxxz, g_x_0_xxxz_xxy, g_x_0_xxxz_xxyz, g_x_0_xxxz_xxz, g_x_0_xxxz_xxzz, g_x_0_xxxz_xyy, g_x_0_xxxz_xyyz, g_x_0_xxxz_xyz, g_x_0_xxxz_xyzz, g_x_0_xxxz_xzz, g_x_0_xxxz_xzzz, g_x_0_xxxz_yyy, g_x_0_xxxz_yyyz, g_x_0_xxxz_yyz, g_x_0_xxxz_yyzz, g_x_0_xxxz_yzz, g_x_0_xxxz_yzzz, g_x_0_xxxz_zzz, g_x_0_xxxz_zzzz, g_x_0_xxxzz_xxx, g_x_0_xxxzz_xxy, g_x_0_xxxzz_xxz, g_x_0_xxxzz_xyy, g_x_0_xxxzz_xyz, g_x_0_xxxzz_xzz, g_x_0_xxxzz_yyy, g_x_0_xxxzz_yyz, g_x_0_xxxzz_yzz, g_x_0_xxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_xxx[k] = -g_x_0_xxxz_xxx[k] * cd_z[k] + g_x_0_xxxz_xxxz[k];

                g_x_0_xxxzz_xxy[k] = -g_x_0_xxxz_xxy[k] * cd_z[k] + g_x_0_xxxz_xxyz[k];

                g_x_0_xxxzz_xxz[k] = -g_x_0_xxxz_xxz[k] * cd_z[k] + g_x_0_xxxz_xxzz[k];

                g_x_0_xxxzz_xyy[k] = -g_x_0_xxxz_xyy[k] * cd_z[k] + g_x_0_xxxz_xyyz[k];

                g_x_0_xxxzz_xyz[k] = -g_x_0_xxxz_xyz[k] * cd_z[k] + g_x_0_xxxz_xyzz[k];

                g_x_0_xxxzz_xzz[k] = -g_x_0_xxxz_xzz[k] * cd_z[k] + g_x_0_xxxz_xzzz[k];

                g_x_0_xxxzz_yyy[k] = -g_x_0_xxxz_yyy[k] * cd_z[k] + g_x_0_xxxz_yyyz[k];

                g_x_0_xxxzz_yyz[k] = -g_x_0_xxxz_yyz[k] * cd_z[k] + g_x_0_xxxz_yyzz[k];

                g_x_0_xxxzz_yzz[k] = -g_x_0_xxxz_yzz[k] * cd_z[k] + g_x_0_xxxz_yzzz[k];

                g_x_0_xxxzz_zzz[k] = -g_x_0_xxxz_zzz[k] * cd_z[k] + g_x_0_xxxz_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 69);

            #pragma omp simd aligned(cd_y, g_x_0_xxyy_xxx, g_x_0_xxyy_xxxy, g_x_0_xxyy_xxy, g_x_0_xxyy_xxyy, g_x_0_xxyy_xxyz, g_x_0_xxyy_xxz, g_x_0_xxyy_xyy, g_x_0_xxyy_xyyy, g_x_0_xxyy_xyyz, g_x_0_xxyy_xyz, g_x_0_xxyy_xyzz, g_x_0_xxyy_xzz, g_x_0_xxyy_yyy, g_x_0_xxyy_yyyy, g_x_0_xxyy_yyyz, g_x_0_xxyy_yyz, g_x_0_xxyy_yyzz, g_x_0_xxyy_yzz, g_x_0_xxyy_yzzz, g_x_0_xxyy_zzz, g_x_0_xxyyy_xxx, g_x_0_xxyyy_xxy, g_x_0_xxyyy_xxz, g_x_0_xxyyy_xyy, g_x_0_xxyyy_xyz, g_x_0_xxyyy_xzz, g_x_0_xxyyy_yyy, g_x_0_xxyyy_yyz, g_x_0_xxyyy_yzz, g_x_0_xxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_xxx[k] = -g_x_0_xxyy_xxx[k] * cd_y[k] + g_x_0_xxyy_xxxy[k];

                g_x_0_xxyyy_xxy[k] = -g_x_0_xxyy_xxy[k] * cd_y[k] + g_x_0_xxyy_xxyy[k];

                g_x_0_xxyyy_xxz[k] = -g_x_0_xxyy_xxz[k] * cd_y[k] + g_x_0_xxyy_xxyz[k];

                g_x_0_xxyyy_xyy[k] = -g_x_0_xxyy_xyy[k] * cd_y[k] + g_x_0_xxyy_xyyy[k];

                g_x_0_xxyyy_xyz[k] = -g_x_0_xxyy_xyz[k] * cd_y[k] + g_x_0_xxyy_xyyz[k];

                g_x_0_xxyyy_xzz[k] = -g_x_0_xxyy_xzz[k] * cd_y[k] + g_x_0_xxyy_xyzz[k];

                g_x_0_xxyyy_yyy[k] = -g_x_0_xxyy_yyy[k] * cd_y[k] + g_x_0_xxyy_yyyy[k];

                g_x_0_xxyyy_yyz[k] = -g_x_0_xxyy_yyz[k] * cd_y[k] + g_x_0_xxyy_yyyz[k];

                g_x_0_xxyyy_yzz[k] = -g_x_0_xxyy_yzz[k] * cd_y[k] + g_x_0_xxyy_yyzz[k];

                g_x_0_xxyyy_zzz[k] = -g_x_0_xxyy_zzz[k] * cd_y[k] + g_x_0_xxyy_yzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 79);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyz_xxx, g_x_0_xxyyz_xxy, g_x_0_xxyyz_xxz, g_x_0_xxyyz_xyy, g_x_0_xxyyz_xyz, g_x_0_xxyyz_xzz, g_x_0_xxyyz_yyy, g_x_0_xxyyz_yyz, g_x_0_xxyyz_yzz, g_x_0_xxyyz_zzz, g_x_0_xxyz_xxx, g_x_0_xxyz_xxxy, g_x_0_xxyz_xxy, g_x_0_xxyz_xxyy, g_x_0_xxyz_xxyz, g_x_0_xxyz_xxz, g_x_0_xxyz_xyy, g_x_0_xxyz_xyyy, g_x_0_xxyz_xyyz, g_x_0_xxyz_xyz, g_x_0_xxyz_xyzz, g_x_0_xxyz_xzz, g_x_0_xxyz_yyy, g_x_0_xxyz_yyyy, g_x_0_xxyz_yyyz, g_x_0_xxyz_yyz, g_x_0_xxyz_yyzz, g_x_0_xxyz_yzz, g_x_0_xxyz_yzzz, g_x_0_xxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_xxx[k] = -g_x_0_xxyz_xxx[k] * cd_y[k] + g_x_0_xxyz_xxxy[k];

                g_x_0_xxyyz_xxy[k] = -g_x_0_xxyz_xxy[k] * cd_y[k] + g_x_0_xxyz_xxyy[k];

                g_x_0_xxyyz_xxz[k] = -g_x_0_xxyz_xxz[k] * cd_y[k] + g_x_0_xxyz_xxyz[k];

                g_x_0_xxyyz_xyy[k] = -g_x_0_xxyz_xyy[k] * cd_y[k] + g_x_0_xxyz_xyyy[k];

                g_x_0_xxyyz_xyz[k] = -g_x_0_xxyz_xyz[k] * cd_y[k] + g_x_0_xxyz_xyyz[k];

                g_x_0_xxyyz_xzz[k] = -g_x_0_xxyz_xzz[k] * cd_y[k] + g_x_0_xxyz_xyzz[k];

                g_x_0_xxyyz_yyy[k] = -g_x_0_xxyz_yyy[k] * cd_y[k] + g_x_0_xxyz_yyyy[k];

                g_x_0_xxyyz_yyz[k] = -g_x_0_xxyz_yyz[k] * cd_y[k] + g_x_0_xxyz_yyyz[k];

                g_x_0_xxyyz_yzz[k] = -g_x_0_xxyz_yzz[k] * cd_y[k] + g_x_0_xxyz_yyzz[k];

                g_x_0_xxyyz_zzz[k] = -g_x_0_xxyz_zzz[k] * cd_y[k] + g_x_0_xxyz_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzz_xxx, g_x_0_xxyzz_xxy, g_x_0_xxyzz_xxz, g_x_0_xxyzz_xyy, g_x_0_xxyzz_xyz, g_x_0_xxyzz_xzz, g_x_0_xxyzz_yyy, g_x_0_xxyzz_yyz, g_x_0_xxyzz_yzz, g_x_0_xxyzz_zzz, g_x_0_xxzz_xxx, g_x_0_xxzz_xxxy, g_x_0_xxzz_xxy, g_x_0_xxzz_xxyy, g_x_0_xxzz_xxyz, g_x_0_xxzz_xxz, g_x_0_xxzz_xyy, g_x_0_xxzz_xyyy, g_x_0_xxzz_xyyz, g_x_0_xxzz_xyz, g_x_0_xxzz_xyzz, g_x_0_xxzz_xzz, g_x_0_xxzz_yyy, g_x_0_xxzz_yyyy, g_x_0_xxzz_yyyz, g_x_0_xxzz_yyz, g_x_0_xxzz_yyzz, g_x_0_xxzz_yzz, g_x_0_xxzz_yzzz, g_x_0_xxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_xxx[k] = -g_x_0_xxzz_xxx[k] * cd_y[k] + g_x_0_xxzz_xxxy[k];

                g_x_0_xxyzz_xxy[k] = -g_x_0_xxzz_xxy[k] * cd_y[k] + g_x_0_xxzz_xxyy[k];

                g_x_0_xxyzz_xxz[k] = -g_x_0_xxzz_xxz[k] * cd_y[k] + g_x_0_xxzz_xxyz[k];

                g_x_0_xxyzz_xyy[k] = -g_x_0_xxzz_xyy[k] * cd_y[k] + g_x_0_xxzz_xyyy[k];

                g_x_0_xxyzz_xyz[k] = -g_x_0_xxzz_xyz[k] * cd_y[k] + g_x_0_xxzz_xyyz[k];

                g_x_0_xxyzz_xzz[k] = -g_x_0_xxzz_xzz[k] * cd_y[k] + g_x_0_xxzz_xyzz[k];

                g_x_0_xxyzz_yyy[k] = -g_x_0_xxzz_yyy[k] * cd_y[k] + g_x_0_xxzz_yyyy[k];

                g_x_0_xxyzz_yyz[k] = -g_x_0_xxzz_yyz[k] * cd_y[k] + g_x_0_xxzz_yyyz[k];

                g_x_0_xxyzz_yzz[k] = -g_x_0_xxzz_yzz[k] * cd_y[k] + g_x_0_xxzz_yyzz[k];

                g_x_0_xxyzz_zzz[k] = -g_x_0_xxzz_zzz[k] * cd_y[k] + g_x_0_xxzz_yzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 99);

            #pragma omp simd aligned(cd_z, g_x_0_xxzz_xxx, g_x_0_xxzz_xxxz, g_x_0_xxzz_xxy, g_x_0_xxzz_xxyz, g_x_0_xxzz_xxz, g_x_0_xxzz_xxzz, g_x_0_xxzz_xyy, g_x_0_xxzz_xyyz, g_x_0_xxzz_xyz, g_x_0_xxzz_xyzz, g_x_0_xxzz_xzz, g_x_0_xxzz_xzzz, g_x_0_xxzz_yyy, g_x_0_xxzz_yyyz, g_x_0_xxzz_yyz, g_x_0_xxzz_yyzz, g_x_0_xxzz_yzz, g_x_0_xxzz_yzzz, g_x_0_xxzz_zzz, g_x_0_xxzz_zzzz, g_x_0_xxzzz_xxx, g_x_0_xxzzz_xxy, g_x_0_xxzzz_xxz, g_x_0_xxzzz_xyy, g_x_0_xxzzz_xyz, g_x_0_xxzzz_xzz, g_x_0_xxzzz_yyy, g_x_0_xxzzz_yyz, g_x_0_xxzzz_yzz, g_x_0_xxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_xxx[k] = -g_x_0_xxzz_xxx[k] * cd_z[k] + g_x_0_xxzz_xxxz[k];

                g_x_0_xxzzz_xxy[k] = -g_x_0_xxzz_xxy[k] * cd_z[k] + g_x_0_xxzz_xxyz[k];

                g_x_0_xxzzz_xxz[k] = -g_x_0_xxzz_xxz[k] * cd_z[k] + g_x_0_xxzz_xxzz[k];

                g_x_0_xxzzz_xyy[k] = -g_x_0_xxzz_xyy[k] * cd_z[k] + g_x_0_xxzz_xyyz[k];

                g_x_0_xxzzz_xyz[k] = -g_x_0_xxzz_xyz[k] * cd_z[k] + g_x_0_xxzz_xyzz[k];

                g_x_0_xxzzz_xzz[k] = -g_x_0_xxzz_xzz[k] * cd_z[k] + g_x_0_xxzz_xzzz[k];

                g_x_0_xxzzz_yyy[k] = -g_x_0_xxzz_yyy[k] * cd_z[k] + g_x_0_xxzz_yyyz[k];

                g_x_0_xxzzz_yyz[k] = -g_x_0_xxzz_yyz[k] * cd_z[k] + g_x_0_xxzz_yyzz[k];

                g_x_0_xxzzz_yzz[k] = -g_x_0_xxzz_yzz[k] * cd_z[k] + g_x_0_xxzz_yzzz[k];

                g_x_0_xxzzz_zzz[k] = -g_x_0_xxzz_zzz[k] * cd_z[k] + g_x_0_xxzz_zzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 109);

            #pragma omp simd aligned(cd_y, g_x_0_xyyy_xxx, g_x_0_xyyy_xxxy, g_x_0_xyyy_xxy, g_x_0_xyyy_xxyy, g_x_0_xyyy_xxyz, g_x_0_xyyy_xxz, g_x_0_xyyy_xyy, g_x_0_xyyy_xyyy, g_x_0_xyyy_xyyz, g_x_0_xyyy_xyz, g_x_0_xyyy_xyzz, g_x_0_xyyy_xzz, g_x_0_xyyy_yyy, g_x_0_xyyy_yyyy, g_x_0_xyyy_yyyz, g_x_0_xyyy_yyz, g_x_0_xyyy_yyzz, g_x_0_xyyy_yzz, g_x_0_xyyy_yzzz, g_x_0_xyyy_zzz, g_x_0_xyyyy_xxx, g_x_0_xyyyy_xxy, g_x_0_xyyyy_xxz, g_x_0_xyyyy_xyy, g_x_0_xyyyy_xyz, g_x_0_xyyyy_xzz, g_x_0_xyyyy_yyy, g_x_0_xyyyy_yyz, g_x_0_xyyyy_yzz, g_x_0_xyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_xxx[k] = -g_x_0_xyyy_xxx[k] * cd_y[k] + g_x_0_xyyy_xxxy[k];

                g_x_0_xyyyy_xxy[k] = -g_x_0_xyyy_xxy[k] * cd_y[k] + g_x_0_xyyy_xxyy[k];

                g_x_0_xyyyy_xxz[k] = -g_x_0_xyyy_xxz[k] * cd_y[k] + g_x_0_xyyy_xxyz[k];

                g_x_0_xyyyy_xyy[k] = -g_x_0_xyyy_xyy[k] * cd_y[k] + g_x_0_xyyy_xyyy[k];

                g_x_0_xyyyy_xyz[k] = -g_x_0_xyyy_xyz[k] * cd_y[k] + g_x_0_xyyy_xyyz[k];

                g_x_0_xyyyy_xzz[k] = -g_x_0_xyyy_xzz[k] * cd_y[k] + g_x_0_xyyy_xyzz[k];

                g_x_0_xyyyy_yyy[k] = -g_x_0_xyyy_yyy[k] * cd_y[k] + g_x_0_xyyy_yyyy[k];

                g_x_0_xyyyy_yyz[k] = -g_x_0_xyyy_yyz[k] * cd_y[k] + g_x_0_xyyy_yyyz[k];

                g_x_0_xyyyy_yzz[k] = -g_x_0_xyyy_yzz[k] * cd_y[k] + g_x_0_xyyy_yyzz[k];

                g_x_0_xyyyy_zzz[k] = -g_x_0_xyyy_zzz[k] * cd_y[k] + g_x_0_xyyy_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyz_xxx, g_x_0_xyyyz_xxy, g_x_0_xyyyz_xxz, g_x_0_xyyyz_xyy, g_x_0_xyyyz_xyz, g_x_0_xyyyz_xzz, g_x_0_xyyyz_yyy, g_x_0_xyyyz_yyz, g_x_0_xyyyz_yzz, g_x_0_xyyyz_zzz, g_x_0_xyyz_xxx, g_x_0_xyyz_xxxy, g_x_0_xyyz_xxy, g_x_0_xyyz_xxyy, g_x_0_xyyz_xxyz, g_x_0_xyyz_xxz, g_x_0_xyyz_xyy, g_x_0_xyyz_xyyy, g_x_0_xyyz_xyyz, g_x_0_xyyz_xyz, g_x_0_xyyz_xyzz, g_x_0_xyyz_xzz, g_x_0_xyyz_yyy, g_x_0_xyyz_yyyy, g_x_0_xyyz_yyyz, g_x_0_xyyz_yyz, g_x_0_xyyz_yyzz, g_x_0_xyyz_yzz, g_x_0_xyyz_yzzz, g_x_0_xyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_xxx[k] = -g_x_0_xyyz_xxx[k] * cd_y[k] + g_x_0_xyyz_xxxy[k];

                g_x_0_xyyyz_xxy[k] = -g_x_0_xyyz_xxy[k] * cd_y[k] + g_x_0_xyyz_xxyy[k];

                g_x_0_xyyyz_xxz[k] = -g_x_0_xyyz_xxz[k] * cd_y[k] + g_x_0_xyyz_xxyz[k];

                g_x_0_xyyyz_xyy[k] = -g_x_0_xyyz_xyy[k] * cd_y[k] + g_x_0_xyyz_xyyy[k];

                g_x_0_xyyyz_xyz[k] = -g_x_0_xyyz_xyz[k] * cd_y[k] + g_x_0_xyyz_xyyz[k];

                g_x_0_xyyyz_xzz[k] = -g_x_0_xyyz_xzz[k] * cd_y[k] + g_x_0_xyyz_xyzz[k];

                g_x_0_xyyyz_yyy[k] = -g_x_0_xyyz_yyy[k] * cd_y[k] + g_x_0_xyyz_yyyy[k];

                g_x_0_xyyyz_yyz[k] = -g_x_0_xyyz_yyz[k] * cd_y[k] + g_x_0_xyyz_yyyz[k];

                g_x_0_xyyyz_yzz[k] = -g_x_0_xyyz_yzz[k] * cd_y[k] + g_x_0_xyyz_yyzz[k];

                g_x_0_xyyyz_zzz[k] = -g_x_0_xyyz_zzz[k] * cd_y[k] + g_x_0_xyyz_yzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 129);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzz_xxx, g_x_0_xyyzz_xxy, g_x_0_xyyzz_xxz, g_x_0_xyyzz_xyy, g_x_0_xyyzz_xyz, g_x_0_xyyzz_xzz, g_x_0_xyyzz_yyy, g_x_0_xyyzz_yyz, g_x_0_xyyzz_yzz, g_x_0_xyyzz_zzz, g_x_0_xyzz_xxx, g_x_0_xyzz_xxxy, g_x_0_xyzz_xxy, g_x_0_xyzz_xxyy, g_x_0_xyzz_xxyz, g_x_0_xyzz_xxz, g_x_0_xyzz_xyy, g_x_0_xyzz_xyyy, g_x_0_xyzz_xyyz, g_x_0_xyzz_xyz, g_x_0_xyzz_xyzz, g_x_0_xyzz_xzz, g_x_0_xyzz_yyy, g_x_0_xyzz_yyyy, g_x_0_xyzz_yyyz, g_x_0_xyzz_yyz, g_x_0_xyzz_yyzz, g_x_0_xyzz_yzz, g_x_0_xyzz_yzzz, g_x_0_xyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_xxx[k] = -g_x_0_xyzz_xxx[k] * cd_y[k] + g_x_0_xyzz_xxxy[k];

                g_x_0_xyyzz_xxy[k] = -g_x_0_xyzz_xxy[k] * cd_y[k] + g_x_0_xyzz_xxyy[k];

                g_x_0_xyyzz_xxz[k] = -g_x_0_xyzz_xxz[k] * cd_y[k] + g_x_0_xyzz_xxyz[k];

                g_x_0_xyyzz_xyy[k] = -g_x_0_xyzz_xyy[k] * cd_y[k] + g_x_0_xyzz_xyyy[k];

                g_x_0_xyyzz_xyz[k] = -g_x_0_xyzz_xyz[k] * cd_y[k] + g_x_0_xyzz_xyyz[k];

                g_x_0_xyyzz_xzz[k] = -g_x_0_xyzz_xzz[k] * cd_y[k] + g_x_0_xyzz_xyzz[k];

                g_x_0_xyyzz_yyy[k] = -g_x_0_xyzz_yyy[k] * cd_y[k] + g_x_0_xyzz_yyyy[k];

                g_x_0_xyyzz_yyz[k] = -g_x_0_xyzz_yyz[k] * cd_y[k] + g_x_0_xyzz_yyyz[k];

                g_x_0_xyyzz_yzz[k] = -g_x_0_xyzz_yzz[k] * cd_y[k] + g_x_0_xyzz_yyzz[k];

                g_x_0_xyyzz_zzz[k] = -g_x_0_xyzz_zzz[k] * cd_y[k] + g_x_0_xyzz_yzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzz_xxx, g_x_0_xyzzz_xxy, g_x_0_xyzzz_xxz, g_x_0_xyzzz_xyy, g_x_0_xyzzz_xyz, g_x_0_xyzzz_xzz, g_x_0_xyzzz_yyy, g_x_0_xyzzz_yyz, g_x_0_xyzzz_yzz, g_x_0_xyzzz_zzz, g_x_0_xzzz_xxx, g_x_0_xzzz_xxxy, g_x_0_xzzz_xxy, g_x_0_xzzz_xxyy, g_x_0_xzzz_xxyz, g_x_0_xzzz_xxz, g_x_0_xzzz_xyy, g_x_0_xzzz_xyyy, g_x_0_xzzz_xyyz, g_x_0_xzzz_xyz, g_x_0_xzzz_xyzz, g_x_0_xzzz_xzz, g_x_0_xzzz_yyy, g_x_0_xzzz_yyyy, g_x_0_xzzz_yyyz, g_x_0_xzzz_yyz, g_x_0_xzzz_yyzz, g_x_0_xzzz_yzz, g_x_0_xzzz_yzzz, g_x_0_xzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_xxx[k] = -g_x_0_xzzz_xxx[k] * cd_y[k] + g_x_0_xzzz_xxxy[k];

                g_x_0_xyzzz_xxy[k] = -g_x_0_xzzz_xxy[k] * cd_y[k] + g_x_0_xzzz_xxyy[k];

                g_x_0_xyzzz_xxz[k] = -g_x_0_xzzz_xxz[k] * cd_y[k] + g_x_0_xzzz_xxyz[k];

                g_x_0_xyzzz_xyy[k] = -g_x_0_xzzz_xyy[k] * cd_y[k] + g_x_0_xzzz_xyyy[k];

                g_x_0_xyzzz_xyz[k] = -g_x_0_xzzz_xyz[k] * cd_y[k] + g_x_0_xzzz_xyyz[k];

                g_x_0_xyzzz_xzz[k] = -g_x_0_xzzz_xzz[k] * cd_y[k] + g_x_0_xzzz_xyzz[k];

                g_x_0_xyzzz_yyy[k] = -g_x_0_xzzz_yyy[k] * cd_y[k] + g_x_0_xzzz_yyyy[k];

                g_x_0_xyzzz_yyz[k] = -g_x_0_xzzz_yyz[k] * cd_y[k] + g_x_0_xzzz_yyyz[k];

                g_x_0_xyzzz_yzz[k] = -g_x_0_xzzz_yzz[k] * cd_y[k] + g_x_0_xzzz_yyzz[k];

                g_x_0_xyzzz_zzz[k] = -g_x_0_xzzz_zzz[k] * cd_y[k] + g_x_0_xzzz_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_z, g_x_0_xzzz_xxx, g_x_0_xzzz_xxxz, g_x_0_xzzz_xxy, g_x_0_xzzz_xxyz, g_x_0_xzzz_xxz, g_x_0_xzzz_xxzz, g_x_0_xzzz_xyy, g_x_0_xzzz_xyyz, g_x_0_xzzz_xyz, g_x_0_xzzz_xyzz, g_x_0_xzzz_xzz, g_x_0_xzzz_xzzz, g_x_0_xzzz_yyy, g_x_0_xzzz_yyyz, g_x_0_xzzz_yyz, g_x_0_xzzz_yyzz, g_x_0_xzzz_yzz, g_x_0_xzzz_yzzz, g_x_0_xzzz_zzz, g_x_0_xzzz_zzzz, g_x_0_xzzzz_xxx, g_x_0_xzzzz_xxy, g_x_0_xzzzz_xxz, g_x_0_xzzzz_xyy, g_x_0_xzzzz_xyz, g_x_0_xzzzz_xzz, g_x_0_xzzzz_yyy, g_x_0_xzzzz_yyz, g_x_0_xzzzz_yzz, g_x_0_xzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_xxx[k] = -g_x_0_xzzz_xxx[k] * cd_z[k] + g_x_0_xzzz_xxxz[k];

                g_x_0_xzzzz_xxy[k] = -g_x_0_xzzz_xxy[k] * cd_z[k] + g_x_0_xzzz_xxyz[k];

                g_x_0_xzzzz_xxz[k] = -g_x_0_xzzz_xxz[k] * cd_z[k] + g_x_0_xzzz_xxzz[k];

                g_x_0_xzzzz_xyy[k] = -g_x_0_xzzz_xyy[k] * cd_z[k] + g_x_0_xzzz_xyyz[k];

                g_x_0_xzzzz_xyz[k] = -g_x_0_xzzz_xyz[k] * cd_z[k] + g_x_0_xzzz_xyzz[k];

                g_x_0_xzzzz_xzz[k] = -g_x_0_xzzz_xzz[k] * cd_z[k] + g_x_0_xzzz_xzzz[k];

                g_x_0_xzzzz_yyy[k] = -g_x_0_xzzz_yyy[k] * cd_z[k] + g_x_0_xzzz_yyyz[k];

                g_x_0_xzzzz_yyz[k] = -g_x_0_xzzz_yyz[k] * cd_z[k] + g_x_0_xzzz_yyzz[k];

                g_x_0_xzzzz_yzz[k] = -g_x_0_xzzz_yzz[k] * cd_z[k] + g_x_0_xzzz_yzzz[k];

                g_x_0_xzzzz_zzz[k] = -g_x_0_xzzz_zzz[k] * cd_z[k] + g_x_0_xzzz_zzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 159);

            #pragma omp simd aligned(cd_y, g_x_0_yyyy_xxx, g_x_0_yyyy_xxxy, g_x_0_yyyy_xxy, g_x_0_yyyy_xxyy, g_x_0_yyyy_xxyz, g_x_0_yyyy_xxz, g_x_0_yyyy_xyy, g_x_0_yyyy_xyyy, g_x_0_yyyy_xyyz, g_x_0_yyyy_xyz, g_x_0_yyyy_xyzz, g_x_0_yyyy_xzz, g_x_0_yyyy_yyy, g_x_0_yyyy_yyyy, g_x_0_yyyy_yyyz, g_x_0_yyyy_yyz, g_x_0_yyyy_yyzz, g_x_0_yyyy_yzz, g_x_0_yyyy_yzzz, g_x_0_yyyy_zzz, g_x_0_yyyyy_xxx, g_x_0_yyyyy_xxy, g_x_0_yyyyy_xxz, g_x_0_yyyyy_xyy, g_x_0_yyyyy_xyz, g_x_0_yyyyy_xzz, g_x_0_yyyyy_yyy, g_x_0_yyyyy_yyz, g_x_0_yyyyy_yzz, g_x_0_yyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_xxx[k] = -g_x_0_yyyy_xxx[k] * cd_y[k] + g_x_0_yyyy_xxxy[k];

                g_x_0_yyyyy_xxy[k] = -g_x_0_yyyy_xxy[k] * cd_y[k] + g_x_0_yyyy_xxyy[k];

                g_x_0_yyyyy_xxz[k] = -g_x_0_yyyy_xxz[k] * cd_y[k] + g_x_0_yyyy_xxyz[k];

                g_x_0_yyyyy_xyy[k] = -g_x_0_yyyy_xyy[k] * cd_y[k] + g_x_0_yyyy_xyyy[k];

                g_x_0_yyyyy_xyz[k] = -g_x_0_yyyy_xyz[k] * cd_y[k] + g_x_0_yyyy_xyyz[k];

                g_x_0_yyyyy_xzz[k] = -g_x_0_yyyy_xzz[k] * cd_y[k] + g_x_0_yyyy_xyzz[k];

                g_x_0_yyyyy_yyy[k] = -g_x_0_yyyy_yyy[k] * cd_y[k] + g_x_0_yyyy_yyyy[k];

                g_x_0_yyyyy_yyz[k] = -g_x_0_yyyy_yyz[k] * cd_y[k] + g_x_0_yyyy_yyyz[k];

                g_x_0_yyyyy_yzz[k] = -g_x_0_yyyy_yzz[k] * cd_y[k] + g_x_0_yyyy_yyzz[k];

                g_x_0_yyyyy_zzz[k] = -g_x_0_yyyy_zzz[k] * cd_y[k] + g_x_0_yyyy_yzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 169);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyz_xxx, g_x_0_yyyyz_xxy, g_x_0_yyyyz_xxz, g_x_0_yyyyz_xyy, g_x_0_yyyyz_xyz, g_x_0_yyyyz_xzz, g_x_0_yyyyz_yyy, g_x_0_yyyyz_yyz, g_x_0_yyyyz_yzz, g_x_0_yyyyz_zzz, g_x_0_yyyz_xxx, g_x_0_yyyz_xxxy, g_x_0_yyyz_xxy, g_x_0_yyyz_xxyy, g_x_0_yyyz_xxyz, g_x_0_yyyz_xxz, g_x_0_yyyz_xyy, g_x_0_yyyz_xyyy, g_x_0_yyyz_xyyz, g_x_0_yyyz_xyz, g_x_0_yyyz_xyzz, g_x_0_yyyz_xzz, g_x_0_yyyz_yyy, g_x_0_yyyz_yyyy, g_x_0_yyyz_yyyz, g_x_0_yyyz_yyz, g_x_0_yyyz_yyzz, g_x_0_yyyz_yzz, g_x_0_yyyz_yzzz, g_x_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_xxx[k] = -g_x_0_yyyz_xxx[k] * cd_y[k] + g_x_0_yyyz_xxxy[k];

                g_x_0_yyyyz_xxy[k] = -g_x_0_yyyz_xxy[k] * cd_y[k] + g_x_0_yyyz_xxyy[k];

                g_x_0_yyyyz_xxz[k] = -g_x_0_yyyz_xxz[k] * cd_y[k] + g_x_0_yyyz_xxyz[k];

                g_x_0_yyyyz_xyy[k] = -g_x_0_yyyz_xyy[k] * cd_y[k] + g_x_0_yyyz_xyyy[k];

                g_x_0_yyyyz_xyz[k] = -g_x_0_yyyz_xyz[k] * cd_y[k] + g_x_0_yyyz_xyyz[k];

                g_x_0_yyyyz_xzz[k] = -g_x_0_yyyz_xzz[k] * cd_y[k] + g_x_0_yyyz_xyzz[k];

                g_x_0_yyyyz_yyy[k] = -g_x_0_yyyz_yyy[k] * cd_y[k] + g_x_0_yyyz_yyyy[k];

                g_x_0_yyyyz_yyz[k] = -g_x_0_yyyz_yyz[k] * cd_y[k] + g_x_0_yyyz_yyyz[k];

                g_x_0_yyyyz_yzz[k] = -g_x_0_yyyz_yzz[k] * cd_y[k] + g_x_0_yyyz_yyzz[k];

                g_x_0_yyyyz_zzz[k] = -g_x_0_yyyz_zzz[k] * cd_y[k] + g_x_0_yyyz_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzz_xxx, g_x_0_yyyzz_xxy, g_x_0_yyyzz_xxz, g_x_0_yyyzz_xyy, g_x_0_yyyzz_xyz, g_x_0_yyyzz_xzz, g_x_0_yyyzz_yyy, g_x_0_yyyzz_yyz, g_x_0_yyyzz_yzz, g_x_0_yyyzz_zzz, g_x_0_yyzz_xxx, g_x_0_yyzz_xxxy, g_x_0_yyzz_xxy, g_x_0_yyzz_xxyy, g_x_0_yyzz_xxyz, g_x_0_yyzz_xxz, g_x_0_yyzz_xyy, g_x_0_yyzz_xyyy, g_x_0_yyzz_xyyz, g_x_0_yyzz_xyz, g_x_0_yyzz_xyzz, g_x_0_yyzz_xzz, g_x_0_yyzz_yyy, g_x_0_yyzz_yyyy, g_x_0_yyzz_yyyz, g_x_0_yyzz_yyz, g_x_0_yyzz_yyzz, g_x_0_yyzz_yzz, g_x_0_yyzz_yzzz, g_x_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_xxx[k] = -g_x_0_yyzz_xxx[k] * cd_y[k] + g_x_0_yyzz_xxxy[k];

                g_x_0_yyyzz_xxy[k] = -g_x_0_yyzz_xxy[k] * cd_y[k] + g_x_0_yyzz_xxyy[k];

                g_x_0_yyyzz_xxz[k] = -g_x_0_yyzz_xxz[k] * cd_y[k] + g_x_0_yyzz_xxyz[k];

                g_x_0_yyyzz_xyy[k] = -g_x_0_yyzz_xyy[k] * cd_y[k] + g_x_0_yyzz_xyyy[k];

                g_x_0_yyyzz_xyz[k] = -g_x_0_yyzz_xyz[k] * cd_y[k] + g_x_0_yyzz_xyyz[k];

                g_x_0_yyyzz_xzz[k] = -g_x_0_yyzz_xzz[k] * cd_y[k] + g_x_0_yyzz_xyzz[k];

                g_x_0_yyyzz_yyy[k] = -g_x_0_yyzz_yyy[k] * cd_y[k] + g_x_0_yyzz_yyyy[k];

                g_x_0_yyyzz_yyz[k] = -g_x_0_yyzz_yyz[k] * cd_y[k] + g_x_0_yyzz_yyyz[k];

                g_x_0_yyyzz_yzz[k] = -g_x_0_yyzz_yzz[k] * cd_y[k] + g_x_0_yyzz_yyzz[k];

                g_x_0_yyyzz_zzz[k] = -g_x_0_yyzz_zzz[k] * cd_y[k] + g_x_0_yyzz_yzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 189);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzz_xxx, g_x_0_yyzzz_xxy, g_x_0_yyzzz_xxz, g_x_0_yyzzz_xyy, g_x_0_yyzzz_xyz, g_x_0_yyzzz_xzz, g_x_0_yyzzz_yyy, g_x_0_yyzzz_yyz, g_x_0_yyzzz_yzz, g_x_0_yyzzz_zzz, g_x_0_yzzz_xxx, g_x_0_yzzz_xxxy, g_x_0_yzzz_xxy, g_x_0_yzzz_xxyy, g_x_0_yzzz_xxyz, g_x_0_yzzz_xxz, g_x_0_yzzz_xyy, g_x_0_yzzz_xyyy, g_x_0_yzzz_xyyz, g_x_0_yzzz_xyz, g_x_0_yzzz_xyzz, g_x_0_yzzz_xzz, g_x_0_yzzz_yyy, g_x_0_yzzz_yyyy, g_x_0_yzzz_yyyz, g_x_0_yzzz_yyz, g_x_0_yzzz_yyzz, g_x_0_yzzz_yzz, g_x_0_yzzz_yzzz, g_x_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_xxx[k] = -g_x_0_yzzz_xxx[k] * cd_y[k] + g_x_0_yzzz_xxxy[k];

                g_x_0_yyzzz_xxy[k] = -g_x_0_yzzz_xxy[k] * cd_y[k] + g_x_0_yzzz_xxyy[k];

                g_x_0_yyzzz_xxz[k] = -g_x_0_yzzz_xxz[k] * cd_y[k] + g_x_0_yzzz_xxyz[k];

                g_x_0_yyzzz_xyy[k] = -g_x_0_yzzz_xyy[k] * cd_y[k] + g_x_0_yzzz_xyyy[k];

                g_x_0_yyzzz_xyz[k] = -g_x_0_yzzz_xyz[k] * cd_y[k] + g_x_0_yzzz_xyyz[k];

                g_x_0_yyzzz_xzz[k] = -g_x_0_yzzz_xzz[k] * cd_y[k] + g_x_0_yzzz_xyzz[k];

                g_x_0_yyzzz_yyy[k] = -g_x_0_yzzz_yyy[k] * cd_y[k] + g_x_0_yzzz_yyyy[k];

                g_x_0_yyzzz_yyz[k] = -g_x_0_yzzz_yyz[k] * cd_y[k] + g_x_0_yzzz_yyyz[k];

                g_x_0_yyzzz_yzz[k] = -g_x_0_yzzz_yzz[k] * cd_y[k] + g_x_0_yzzz_yyzz[k];

                g_x_0_yyzzz_zzz[k] = -g_x_0_yzzz_zzz[k] * cd_y[k] + g_x_0_yzzz_yzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 199);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzz_xxx, g_x_0_yzzzz_xxy, g_x_0_yzzzz_xxz, g_x_0_yzzzz_xyy, g_x_0_yzzzz_xyz, g_x_0_yzzzz_xzz, g_x_0_yzzzz_yyy, g_x_0_yzzzz_yyz, g_x_0_yzzzz_yzz, g_x_0_yzzzz_zzz, g_x_0_zzzz_xxx, g_x_0_zzzz_xxxy, g_x_0_zzzz_xxy, g_x_0_zzzz_xxyy, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxz, g_x_0_zzzz_xyy, g_x_0_zzzz_xyyy, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xzz, g_x_0_zzzz_yyy, g_x_0_zzzz_yyyy, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_xxx[k] = -g_x_0_zzzz_xxx[k] * cd_y[k] + g_x_0_zzzz_xxxy[k];

                g_x_0_yzzzz_xxy[k] = -g_x_0_zzzz_xxy[k] * cd_y[k] + g_x_0_zzzz_xxyy[k];

                g_x_0_yzzzz_xxz[k] = -g_x_0_zzzz_xxz[k] * cd_y[k] + g_x_0_zzzz_xxyz[k];

                g_x_0_yzzzz_xyy[k] = -g_x_0_zzzz_xyy[k] * cd_y[k] + g_x_0_zzzz_xyyy[k];

                g_x_0_yzzzz_xyz[k] = -g_x_0_zzzz_xyz[k] * cd_y[k] + g_x_0_zzzz_xyyz[k];

                g_x_0_yzzzz_xzz[k] = -g_x_0_zzzz_xzz[k] * cd_y[k] + g_x_0_zzzz_xyzz[k];

                g_x_0_yzzzz_yyy[k] = -g_x_0_zzzz_yyy[k] * cd_y[k] + g_x_0_zzzz_yyyy[k];

                g_x_0_yzzzz_yyz[k] = -g_x_0_zzzz_yyz[k] * cd_y[k] + g_x_0_zzzz_yyyz[k];

                g_x_0_yzzzz_yzz[k] = -g_x_0_zzzz_yzz[k] * cd_y[k] + g_x_0_zzzz_yyzz[k];

                g_x_0_yzzzz_zzz[k] = -g_x_0_zzzz_zzz[k] * cd_y[k] + g_x_0_zzzz_yzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_z, g_x_0_zzzz_xxx, g_x_0_zzzz_xxxz, g_x_0_zzzz_xxy, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxz, g_x_0_zzzz_xxzz, g_x_0_zzzz_xyy, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xzz, g_x_0_zzzz_xzzz, g_x_0_zzzz_yyy, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_zzz, g_x_0_zzzz_zzzz, g_x_0_zzzzz_xxx, g_x_0_zzzzz_xxy, g_x_0_zzzzz_xxz, g_x_0_zzzzz_xyy, g_x_0_zzzzz_xyz, g_x_0_zzzzz_xzz, g_x_0_zzzzz_yyy, g_x_0_zzzzz_yyz, g_x_0_zzzzz_yzz, g_x_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_xxx[k] = -g_x_0_zzzz_xxx[k] * cd_z[k] + g_x_0_zzzz_xxxz[k];

                g_x_0_zzzzz_xxy[k] = -g_x_0_zzzz_xxy[k] * cd_z[k] + g_x_0_zzzz_xxyz[k];

                g_x_0_zzzzz_xxz[k] = -g_x_0_zzzz_xxz[k] * cd_z[k] + g_x_0_zzzz_xxzz[k];

                g_x_0_zzzzz_xyy[k] = -g_x_0_zzzz_xyy[k] * cd_z[k] + g_x_0_zzzz_xyyz[k];

                g_x_0_zzzzz_xyz[k] = -g_x_0_zzzz_xyz[k] * cd_z[k] + g_x_0_zzzz_xyzz[k];

                g_x_0_zzzzz_xzz[k] = -g_x_0_zzzz_xzz[k] * cd_z[k] + g_x_0_zzzz_xzzz[k];

                g_x_0_zzzzz_yyy[k] = -g_x_0_zzzz_yyy[k] * cd_z[k] + g_x_0_zzzz_yyyz[k];

                g_x_0_zzzzz_yyz[k] = -g_x_0_zzzz_yyz[k] * cd_z[k] + g_x_0_zzzz_yyzz[k];

                g_x_0_zzzzz_yzz[k] = -g_x_0_zzzz_yzz[k] * cd_z[k] + g_x_0_zzzz_yzzz[k];

                g_x_0_zzzzz_zzz[k] = -g_x_0_zzzz_zzz[k] * cd_z[k] + g_x_0_zzzz_zzzz[k];
            }
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 5);

            auto g_y_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 6);

            auto g_y_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 7);

            auto g_y_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 8);

            auto g_y_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_y_0_xxxx_xxx, g_y_0_xxxx_xxxx, g_y_0_xxxx_xxxy, g_y_0_xxxx_xxxz, g_y_0_xxxx_xxy, g_y_0_xxxx_xxyy, g_y_0_xxxx_xxyz, g_y_0_xxxx_xxz, g_y_0_xxxx_xxzz, g_y_0_xxxx_xyy, g_y_0_xxxx_xyyy, g_y_0_xxxx_xyyz, g_y_0_xxxx_xyz, g_y_0_xxxx_xyzz, g_y_0_xxxx_xzz, g_y_0_xxxx_xzzz, g_y_0_xxxx_yyy, g_y_0_xxxx_yyz, g_y_0_xxxx_yzz, g_y_0_xxxx_zzz, g_y_0_xxxxx_xxx, g_y_0_xxxxx_xxy, g_y_0_xxxxx_xxz, g_y_0_xxxxx_xyy, g_y_0_xxxxx_xyz, g_y_0_xxxxx_xzz, g_y_0_xxxxx_yyy, g_y_0_xxxxx_yyz, g_y_0_xxxxx_yzz, g_y_0_xxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_xxx[k] = -g_y_0_xxxx_xxx[k] * cd_x[k] + g_y_0_xxxx_xxxx[k];

                g_y_0_xxxxx_xxy[k] = -g_y_0_xxxx_xxy[k] * cd_x[k] + g_y_0_xxxx_xxxy[k];

                g_y_0_xxxxx_xxz[k] = -g_y_0_xxxx_xxz[k] * cd_x[k] + g_y_0_xxxx_xxxz[k];

                g_y_0_xxxxx_xyy[k] = -g_y_0_xxxx_xyy[k] * cd_x[k] + g_y_0_xxxx_xxyy[k];

                g_y_0_xxxxx_xyz[k] = -g_y_0_xxxx_xyz[k] * cd_x[k] + g_y_0_xxxx_xxyz[k];

                g_y_0_xxxxx_xzz[k] = -g_y_0_xxxx_xzz[k] * cd_x[k] + g_y_0_xxxx_xxzz[k];

                g_y_0_xxxxx_yyy[k] = -g_y_0_xxxx_yyy[k] * cd_x[k] + g_y_0_xxxx_xyyy[k];

                g_y_0_xxxxx_yyz[k] = -g_y_0_xxxx_yyz[k] * cd_x[k] + g_y_0_xxxx_xyyz[k];

                g_y_0_xxxxx_yzz[k] = -g_y_0_xxxx_yzz[k] * cd_x[k] + g_y_0_xxxx_xyzz[k];

                g_y_0_xxxxx_zzz[k] = -g_y_0_xxxx_zzz[k] * cd_x[k] + g_y_0_xxxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 10);

            auto g_y_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 11);

            auto g_y_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 12);

            auto g_y_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 13);

            auto g_y_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 14);

            auto g_y_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 15);

            auto g_y_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 16);

            auto g_y_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 17);

            auto g_y_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 18);

            auto g_y_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxy_xxx, g_y_0_xxxxy_xxy, g_y_0_xxxxy_xxz, g_y_0_xxxxy_xyy, g_y_0_xxxxy_xyz, g_y_0_xxxxy_xzz, g_y_0_xxxxy_yyy, g_y_0_xxxxy_yyz, g_y_0_xxxxy_yzz, g_y_0_xxxxy_zzz, g_y_0_xxxy_xxx, g_y_0_xxxy_xxxx, g_y_0_xxxy_xxxy, g_y_0_xxxy_xxxz, g_y_0_xxxy_xxy, g_y_0_xxxy_xxyy, g_y_0_xxxy_xxyz, g_y_0_xxxy_xxz, g_y_0_xxxy_xxzz, g_y_0_xxxy_xyy, g_y_0_xxxy_xyyy, g_y_0_xxxy_xyyz, g_y_0_xxxy_xyz, g_y_0_xxxy_xyzz, g_y_0_xxxy_xzz, g_y_0_xxxy_xzzz, g_y_0_xxxy_yyy, g_y_0_xxxy_yyz, g_y_0_xxxy_yzz, g_y_0_xxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_xxx[k] = -g_y_0_xxxy_xxx[k] * cd_x[k] + g_y_0_xxxy_xxxx[k];

                g_y_0_xxxxy_xxy[k] = -g_y_0_xxxy_xxy[k] * cd_x[k] + g_y_0_xxxy_xxxy[k];

                g_y_0_xxxxy_xxz[k] = -g_y_0_xxxy_xxz[k] * cd_x[k] + g_y_0_xxxy_xxxz[k];

                g_y_0_xxxxy_xyy[k] = -g_y_0_xxxy_xyy[k] * cd_x[k] + g_y_0_xxxy_xxyy[k];

                g_y_0_xxxxy_xyz[k] = -g_y_0_xxxy_xyz[k] * cd_x[k] + g_y_0_xxxy_xxyz[k];

                g_y_0_xxxxy_xzz[k] = -g_y_0_xxxy_xzz[k] * cd_x[k] + g_y_0_xxxy_xxzz[k];

                g_y_0_xxxxy_yyy[k] = -g_y_0_xxxy_yyy[k] * cd_x[k] + g_y_0_xxxy_xyyy[k];

                g_y_0_xxxxy_yyz[k] = -g_y_0_xxxy_yyz[k] * cd_x[k] + g_y_0_xxxy_xyyz[k];

                g_y_0_xxxxy_yzz[k] = -g_y_0_xxxy_yzz[k] * cd_x[k] + g_y_0_xxxy_xyzz[k];

                g_y_0_xxxxy_zzz[k] = -g_y_0_xxxy_zzz[k] * cd_x[k] + g_y_0_xxxy_xzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 20);

            auto g_y_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 21);

            auto g_y_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 22);

            auto g_y_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 23);

            auto g_y_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 24);

            auto g_y_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 25);

            auto g_y_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 26);

            auto g_y_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 27);

            auto g_y_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 28);

            auto g_y_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxz_xxx, g_y_0_xxxxz_xxy, g_y_0_xxxxz_xxz, g_y_0_xxxxz_xyy, g_y_0_xxxxz_xyz, g_y_0_xxxxz_xzz, g_y_0_xxxxz_yyy, g_y_0_xxxxz_yyz, g_y_0_xxxxz_yzz, g_y_0_xxxxz_zzz, g_y_0_xxxz_xxx, g_y_0_xxxz_xxxx, g_y_0_xxxz_xxxy, g_y_0_xxxz_xxxz, g_y_0_xxxz_xxy, g_y_0_xxxz_xxyy, g_y_0_xxxz_xxyz, g_y_0_xxxz_xxz, g_y_0_xxxz_xxzz, g_y_0_xxxz_xyy, g_y_0_xxxz_xyyy, g_y_0_xxxz_xyyz, g_y_0_xxxz_xyz, g_y_0_xxxz_xyzz, g_y_0_xxxz_xzz, g_y_0_xxxz_xzzz, g_y_0_xxxz_yyy, g_y_0_xxxz_yyz, g_y_0_xxxz_yzz, g_y_0_xxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_xxx[k] = -g_y_0_xxxz_xxx[k] * cd_x[k] + g_y_0_xxxz_xxxx[k];

                g_y_0_xxxxz_xxy[k] = -g_y_0_xxxz_xxy[k] * cd_x[k] + g_y_0_xxxz_xxxy[k];

                g_y_0_xxxxz_xxz[k] = -g_y_0_xxxz_xxz[k] * cd_x[k] + g_y_0_xxxz_xxxz[k];

                g_y_0_xxxxz_xyy[k] = -g_y_0_xxxz_xyy[k] * cd_x[k] + g_y_0_xxxz_xxyy[k];

                g_y_0_xxxxz_xyz[k] = -g_y_0_xxxz_xyz[k] * cd_x[k] + g_y_0_xxxz_xxyz[k];

                g_y_0_xxxxz_xzz[k] = -g_y_0_xxxz_xzz[k] * cd_x[k] + g_y_0_xxxz_xxzz[k];

                g_y_0_xxxxz_yyy[k] = -g_y_0_xxxz_yyy[k] * cd_x[k] + g_y_0_xxxz_xyyy[k];

                g_y_0_xxxxz_yyz[k] = -g_y_0_xxxz_yyz[k] * cd_x[k] + g_y_0_xxxz_xyyz[k];

                g_y_0_xxxxz_yzz[k] = -g_y_0_xxxz_yzz[k] * cd_x[k] + g_y_0_xxxz_xyzz[k];

                g_y_0_xxxxz_zzz[k] = -g_y_0_xxxz_zzz[k] * cd_x[k] + g_y_0_xxxz_xzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 30);

            auto g_y_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 31);

            auto g_y_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 32);

            auto g_y_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 33);

            auto g_y_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 34);

            auto g_y_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 35);

            auto g_y_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 36);

            auto g_y_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 37);

            auto g_y_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 38);

            auto g_y_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 39);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyy_xxx, g_y_0_xxxyy_xxy, g_y_0_xxxyy_xxz, g_y_0_xxxyy_xyy, g_y_0_xxxyy_xyz, g_y_0_xxxyy_xzz, g_y_0_xxxyy_yyy, g_y_0_xxxyy_yyz, g_y_0_xxxyy_yzz, g_y_0_xxxyy_zzz, g_y_0_xxyy_xxx, g_y_0_xxyy_xxxx, g_y_0_xxyy_xxxy, g_y_0_xxyy_xxxz, g_y_0_xxyy_xxy, g_y_0_xxyy_xxyy, g_y_0_xxyy_xxyz, g_y_0_xxyy_xxz, g_y_0_xxyy_xxzz, g_y_0_xxyy_xyy, g_y_0_xxyy_xyyy, g_y_0_xxyy_xyyz, g_y_0_xxyy_xyz, g_y_0_xxyy_xyzz, g_y_0_xxyy_xzz, g_y_0_xxyy_xzzz, g_y_0_xxyy_yyy, g_y_0_xxyy_yyz, g_y_0_xxyy_yzz, g_y_0_xxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_xxx[k] = -g_y_0_xxyy_xxx[k] * cd_x[k] + g_y_0_xxyy_xxxx[k];

                g_y_0_xxxyy_xxy[k] = -g_y_0_xxyy_xxy[k] * cd_x[k] + g_y_0_xxyy_xxxy[k];

                g_y_0_xxxyy_xxz[k] = -g_y_0_xxyy_xxz[k] * cd_x[k] + g_y_0_xxyy_xxxz[k];

                g_y_0_xxxyy_xyy[k] = -g_y_0_xxyy_xyy[k] * cd_x[k] + g_y_0_xxyy_xxyy[k];

                g_y_0_xxxyy_xyz[k] = -g_y_0_xxyy_xyz[k] * cd_x[k] + g_y_0_xxyy_xxyz[k];

                g_y_0_xxxyy_xzz[k] = -g_y_0_xxyy_xzz[k] * cd_x[k] + g_y_0_xxyy_xxzz[k];

                g_y_0_xxxyy_yyy[k] = -g_y_0_xxyy_yyy[k] * cd_x[k] + g_y_0_xxyy_xyyy[k];

                g_y_0_xxxyy_yyz[k] = -g_y_0_xxyy_yyz[k] * cd_x[k] + g_y_0_xxyy_xyyz[k];

                g_y_0_xxxyy_yzz[k] = -g_y_0_xxyy_yzz[k] * cd_x[k] + g_y_0_xxyy_xyzz[k];

                g_y_0_xxxyy_zzz[k] = -g_y_0_xxyy_zzz[k] * cd_x[k] + g_y_0_xxyy_xzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 40);

            auto g_y_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 41);

            auto g_y_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 42);

            auto g_y_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 43);

            auto g_y_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 44);

            auto g_y_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 45);

            auto g_y_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 46);

            auto g_y_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 47);

            auto g_y_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 48);

            auto g_y_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 49);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyz_xxx, g_y_0_xxxyz_xxy, g_y_0_xxxyz_xxz, g_y_0_xxxyz_xyy, g_y_0_xxxyz_xyz, g_y_0_xxxyz_xzz, g_y_0_xxxyz_yyy, g_y_0_xxxyz_yyz, g_y_0_xxxyz_yzz, g_y_0_xxxyz_zzz, g_y_0_xxyz_xxx, g_y_0_xxyz_xxxx, g_y_0_xxyz_xxxy, g_y_0_xxyz_xxxz, g_y_0_xxyz_xxy, g_y_0_xxyz_xxyy, g_y_0_xxyz_xxyz, g_y_0_xxyz_xxz, g_y_0_xxyz_xxzz, g_y_0_xxyz_xyy, g_y_0_xxyz_xyyy, g_y_0_xxyz_xyyz, g_y_0_xxyz_xyz, g_y_0_xxyz_xyzz, g_y_0_xxyz_xzz, g_y_0_xxyz_xzzz, g_y_0_xxyz_yyy, g_y_0_xxyz_yyz, g_y_0_xxyz_yzz, g_y_0_xxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_xxx[k] = -g_y_0_xxyz_xxx[k] * cd_x[k] + g_y_0_xxyz_xxxx[k];

                g_y_0_xxxyz_xxy[k] = -g_y_0_xxyz_xxy[k] * cd_x[k] + g_y_0_xxyz_xxxy[k];

                g_y_0_xxxyz_xxz[k] = -g_y_0_xxyz_xxz[k] * cd_x[k] + g_y_0_xxyz_xxxz[k];

                g_y_0_xxxyz_xyy[k] = -g_y_0_xxyz_xyy[k] * cd_x[k] + g_y_0_xxyz_xxyy[k];

                g_y_0_xxxyz_xyz[k] = -g_y_0_xxyz_xyz[k] * cd_x[k] + g_y_0_xxyz_xxyz[k];

                g_y_0_xxxyz_xzz[k] = -g_y_0_xxyz_xzz[k] * cd_x[k] + g_y_0_xxyz_xxzz[k];

                g_y_0_xxxyz_yyy[k] = -g_y_0_xxyz_yyy[k] * cd_x[k] + g_y_0_xxyz_xyyy[k];

                g_y_0_xxxyz_yyz[k] = -g_y_0_xxyz_yyz[k] * cd_x[k] + g_y_0_xxyz_xyyz[k];

                g_y_0_xxxyz_yzz[k] = -g_y_0_xxyz_yzz[k] * cd_x[k] + g_y_0_xxyz_xyzz[k];

                g_y_0_xxxyz_zzz[k] = -g_y_0_xxyz_zzz[k] * cd_x[k] + g_y_0_xxyz_xzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 50);

            auto g_y_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 51);

            auto g_y_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 52);

            auto g_y_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 53);

            auto g_y_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 54);

            auto g_y_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 55);

            auto g_y_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 56);

            auto g_y_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 57);

            auto g_y_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 58);

            auto g_y_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzz_xxx, g_y_0_xxxzz_xxy, g_y_0_xxxzz_xxz, g_y_0_xxxzz_xyy, g_y_0_xxxzz_xyz, g_y_0_xxxzz_xzz, g_y_0_xxxzz_yyy, g_y_0_xxxzz_yyz, g_y_0_xxxzz_yzz, g_y_0_xxxzz_zzz, g_y_0_xxzz_xxx, g_y_0_xxzz_xxxx, g_y_0_xxzz_xxxy, g_y_0_xxzz_xxxz, g_y_0_xxzz_xxy, g_y_0_xxzz_xxyy, g_y_0_xxzz_xxyz, g_y_0_xxzz_xxz, g_y_0_xxzz_xxzz, g_y_0_xxzz_xyy, g_y_0_xxzz_xyyy, g_y_0_xxzz_xyyz, g_y_0_xxzz_xyz, g_y_0_xxzz_xyzz, g_y_0_xxzz_xzz, g_y_0_xxzz_xzzz, g_y_0_xxzz_yyy, g_y_0_xxzz_yyz, g_y_0_xxzz_yzz, g_y_0_xxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_xxx[k] = -g_y_0_xxzz_xxx[k] * cd_x[k] + g_y_0_xxzz_xxxx[k];

                g_y_0_xxxzz_xxy[k] = -g_y_0_xxzz_xxy[k] * cd_x[k] + g_y_0_xxzz_xxxy[k];

                g_y_0_xxxzz_xxz[k] = -g_y_0_xxzz_xxz[k] * cd_x[k] + g_y_0_xxzz_xxxz[k];

                g_y_0_xxxzz_xyy[k] = -g_y_0_xxzz_xyy[k] * cd_x[k] + g_y_0_xxzz_xxyy[k];

                g_y_0_xxxzz_xyz[k] = -g_y_0_xxzz_xyz[k] * cd_x[k] + g_y_0_xxzz_xxyz[k];

                g_y_0_xxxzz_xzz[k] = -g_y_0_xxzz_xzz[k] * cd_x[k] + g_y_0_xxzz_xxzz[k];

                g_y_0_xxxzz_yyy[k] = -g_y_0_xxzz_yyy[k] * cd_x[k] + g_y_0_xxzz_xyyy[k];

                g_y_0_xxxzz_yyz[k] = -g_y_0_xxzz_yyz[k] * cd_x[k] + g_y_0_xxzz_xyyz[k];

                g_y_0_xxxzz_yzz[k] = -g_y_0_xxzz_yzz[k] * cd_x[k] + g_y_0_xxzz_xyzz[k];

                g_y_0_xxxzz_zzz[k] = -g_y_0_xxzz_zzz[k] * cd_x[k] + g_y_0_xxzz_xzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 60);

            auto g_y_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 61);

            auto g_y_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 62);

            auto g_y_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 63);

            auto g_y_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 64);

            auto g_y_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 65);

            auto g_y_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 66);

            auto g_y_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 67);

            auto g_y_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 68);

            auto g_y_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 69);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyy_xxx, g_y_0_xxyyy_xxy, g_y_0_xxyyy_xxz, g_y_0_xxyyy_xyy, g_y_0_xxyyy_xyz, g_y_0_xxyyy_xzz, g_y_0_xxyyy_yyy, g_y_0_xxyyy_yyz, g_y_0_xxyyy_yzz, g_y_0_xxyyy_zzz, g_y_0_xyyy_xxx, g_y_0_xyyy_xxxx, g_y_0_xyyy_xxxy, g_y_0_xyyy_xxxz, g_y_0_xyyy_xxy, g_y_0_xyyy_xxyy, g_y_0_xyyy_xxyz, g_y_0_xyyy_xxz, g_y_0_xyyy_xxzz, g_y_0_xyyy_xyy, g_y_0_xyyy_xyyy, g_y_0_xyyy_xyyz, g_y_0_xyyy_xyz, g_y_0_xyyy_xyzz, g_y_0_xyyy_xzz, g_y_0_xyyy_xzzz, g_y_0_xyyy_yyy, g_y_0_xyyy_yyz, g_y_0_xyyy_yzz, g_y_0_xyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_xxx[k] = -g_y_0_xyyy_xxx[k] * cd_x[k] + g_y_0_xyyy_xxxx[k];

                g_y_0_xxyyy_xxy[k] = -g_y_0_xyyy_xxy[k] * cd_x[k] + g_y_0_xyyy_xxxy[k];

                g_y_0_xxyyy_xxz[k] = -g_y_0_xyyy_xxz[k] * cd_x[k] + g_y_0_xyyy_xxxz[k];

                g_y_0_xxyyy_xyy[k] = -g_y_0_xyyy_xyy[k] * cd_x[k] + g_y_0_xyyy_xxyy[k];

                g_y_0_xxyyy_xyz[k] = -g_y_0_xyyy_xyz[k] * cd_x[k] + g_y_0_xyyy_xxyz[k];

                g_y_0_xxyyy_xzz[k] = -g_y_0_xyyy_xzz[k] * cd_x[k] + g_y_0_xyyy_xxzz[k];

                g_y_0_xxyyy_yyy[k] = -g_y_0_xyyy_yyy[k] * cd_x[k] + g_y_0_xyyy_xyyy[k];

                g_y_0_xxyyy_yyz[k] = -g_y_0_xyyy_yyz[k] * cd_x[k] + g_y_0_xyyy_xyyz[k];

                g_y_0_xxyyy_yzz[k] = -g_y_0_xyyy_yzz[k] * cd_x[k] + g_y_0_xyyy_xyzz[k];

                g_y_0_xxyyy_zzz[k] = -g_y_0_xyyy_zzz[k] * cd_x[k] + g_y_0_xyyy_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 70);

            auto g_y_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 71);

            auto g_y_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 72);

            auto g_y_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 73);

            auto g_y_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 74);

            auto g_y_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 75);

            auto g_y_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 76);

            auto g_y_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 77);

            auto g_y_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 78);

            auto g_y_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 79);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyz_xxx, g_y_0_xxyyz_xxy, g_y_0_xxyyz_xxz, g_y_0_xxyyz_xyy, g_y_0_xxyyz_xyz, g_y_0_xxyyz_xzz, g_y_0_xxyyz_yyy, g_y_0_xxyyz_yyz, g_y_0_xxyyz_yzz, g_y_0_xxyyz_zzz, g_y_0_xyyz_xxx, g_y_0_xyyz_xxxx, g_y_0_xyyz_xxxy, g_y_0_xyyz_xxxz, g_y_0_xyyz_xxy, g_y_0_xyyz_xxyy, g_y_0_xyyz_xxyz, g_y_0_xyyz_xxz, g_y_0_xyyz_xxzz, g_y_0_xyyz_xyy, g_y_0_xyyz_xyyy, g_y_0_xyyz_xyyz, g_y_0_xyyz_xyz, g_y_0_xyyz_xyzz, g_y_0_xyyz_xzz, g_y_0_xyyz_xzzz, g_y_0_xyyz_yyy, g_y_0_xyyz_yyz, g_y_0_xyyz_yzz, g_y_0_xyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_xxx[k] = -g_y_0_xyyz_xxx[k] * cd_x[k] + g_y_0_xyyz_xxxx[k];

                g_y_0_xxyyz_xxy[k] = -g_y_0_xyyz_xxy[k] * cd_x[k] + g_y_0_xyyz_xxxy[k];

                g_y_0_xxyyz_xxz[k] = -g_y_0_xyyz_xxz[k] * cd_x[k] + g_y_0_xyyz_xxxz[k];

                g_y_0_xxyyz_xyy[k] = -g_y_0_xyyz_xyy[k] * cd_x[k] + g_y_0_xyyz_xxyy[k];

                g_y_0_xxyyz_xyz[k] = -g_y_0_xyyz_xyz[k] * cd_x[k] + g_y_0_xyyz_xxyz[k];

                g_y_0_xxyyz_xzz[k] = -g_y_0_xyyz_xzz[k] * cd_x[k] + g_y_0_xyyz_xxzz[k];

                g_y_0_xxyyz_yyy[k] = -g_y_0_xyyz_yyy[k] * cd_x[k] + g_y_0_xyyz_xyyy[k];

                g_y_0_xxyyz_yyz[k] = -g_y_0_xyyz_yyz[k] * cd_x[k] + g_y_0_xyyz_xyyz[k];

                g_y_0_xxyyz_yzz[k] = -g_y_0_xyyz_yzz[k] * cd_x[k] + g_y_0_xyyz_xyzz[k];

                g_y_0_xxyyz_zzz[k] = -g_y_0_xyyz_zzz[k] * cd_x[k] + g_y_0_xyyz_xzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 80);

            auto g_y_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 81);

            auto g_y_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 82);

            auto g_y_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 83);

            auto g_y_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 84);

            auto g_y_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 85);

            auto g_y_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 86);

            auto g_y_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 87);

            auto g_y_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 88);

            auto g_y_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzz_xxx, g_y_0_xxyzz_xxy, g_y_0_xxyzz_xxz, g_y_0_xxyzz_xyy, g_y_0_xxyzz_xyz, g_y_0_xxyzz_xzz, g_y_0_xxyzz_yyy, g_y_0_xxyzz_yyz, g_y_0_xxyzz_yzz, g_y_0_xxyzz_zzz, g_y_0_xyzz_xxx, g_y_0_xyzz_xxxx, g_y_0_xyzz_xxxy, g_y_0_xyzz_xxxz, g_y_0_xyzz_xxy, g_y_0_xyzz_xxyy, g_y_0_xyzz_xxyz, g_y_0_xyzz_xxz, g_y_0_xyzz_xxzz, g_y_0_xyzz_xyy, g_y_0_xyzz_xyyy, g_y_0_xyzz_xyyz, g_y_0_xyzz_xyz, g_y_0_xyzz_xyzz, g_y_0_xyzz_xzz, g_y_0_xyzz_xzzz, g_y_0_xyzz_yyy, g_y_0_xyzz_yyz, g_y_0_xyzz_yzz, g_y_0_xyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_xxx[k] = -g_y_0_xyzz_xxx[k] * cd_x[k] + g_y_0_xyzz_xxxx[k];

                g_y_0_xxyzz_xxy[k] = -g_y_0_xyzz_xxy[k] * cd_x[k] + g_y_0_xyzz_xxxy[k];

                g_y_0_xxyzz_xxz[k] = -g_y_0_xyzz_xxz[k] * cd_x[k] + g_y_0_xyzz_xxxz[k];

                g_y_0_xxyzz_xyy[k] = -g_y_0_xyzz_xyy[k] * cd_x[k] + g_y_0_xyzz_xxyy[k];

                g_y_0_xxyzz_xyz[k] = -g_y_0_xyzz_xyz[k] * cd_x[k] + g_y_0_xyzz_xxyz[k];

                g_y_0_xxyzz_xzz[k] = -g_y_0_xyzz_xzz[k] * cd_x[k] + g_y_0_xyzz_xxzz[k];

                g_y_0_xxyzz_yyy[k] = -g_y_0_xyzz_yyy[k] * cd_x[k] + g_y_0_xyzz_xyyy[k];

                g_y_0_xxyzz_yyz[k] = -g_y_0_xyzz_yyz[k] * cd_x[k] + g_y_0_xyzz_xyyz[k];

                g_y_0_xxyzz_yzz[k] = -g_y_0_xyzz_yzz[k] * cd_x[k] + g_y_0_xyzz_xyzz[k];

                g_y_0_xxyzz_zzz[k] = -g_y_0_xyzz_zzz[k] * cd_x[k] + g_y_0_xyzz_xzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 90);

            auto g_y_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 91);

            auto g_y_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 92);

            auto g_y_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 93);

            auto g_y_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 94);

            auto g_y_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 95);

            auto g_y_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 96);

            auto g_y_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 97);

            auto g_y_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 98);

            auto g_y_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 99);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzz_xxx, g_y_0_xxzzz_xxy, g_y_0_xxzzz_xxz, g_y_0_xxzzz_xyy, g_y_0_xxzzz_xyz, g_y_0_xxzzz_xzz, g_y_0_xxzzz_yyy, g_y_0_xxzzz_yyz, g_y_0_xxzzz_yzz, g_y_0_xxzzz_zzz, g_y_0_xzzz_xxx, g_y_0_xzzz_xxxx, g_y_0_xzzz_xxxy, g_y_0_xzzz_xxxz, g_y_0_xzzz_xxy, g_y_0_xzzz_xxyy, g_y_0_xzzz_xxyz, g_y_0_xzzz_xxz, g_y_0_xzzz_xxzz, g_y_0_xzzz_xyy, g_y_0_xzzz_xyyy, g_y_0_xzzz_xyyz, g_y_0_xzzz_xyz, g_y_0_xzzz_xyzz, g_y_0_xzzz_xzz, g_y_0_xzzz_xzzz, g_y_0_xzzz_yyy, g_y_0_xzzz_yyz, g_y_0_xzzz_yzz, g_y_0_xzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_xxx[k] = -g_y_0_xzzz_xxx[k] * cd_x[k] + g_y_0_xzzz_xxxx[k];

                g_y_0_xxzzz_xxy[k] = -g_y_0_xzzz_xxy[k] * cd_x[k] + g_y_0_xzzz_xxxy[k];

                g_y_0_xxzzz_xxz[k] = -g_y_0_xzzz_xxz[k] * cd_x[k] + g_y_0_xzzz_xxxz[k];

                g_y_0_xxzzz_xyy[k] = -g_y_0_xzzz_xyy[k] * cd_x[k] + g_y_0_xzzz_xxyy[k];

                g_y_0_xxzzz_xyz[k] = -g_y_0_xzzz_xyz[k] * cd_x[k] + g_y_0_xzzz_xxyz[k];

                g_y_0_xxzzz_xzz[k] = -g_y_0_xzzz_xzz[k] * cd_x[k] + g_y_0_xzzz_xxzz[k];

                g_y_0_xxzzz_yyy[k] = -g_y_0_xzzz_yyy[k] * cd_x[k] + g_y_0_xzzz_xyyy[k];

                g_y_0_xxzzz_yyz[k] = -g_y_0_xzzz_yyz[k] * cd_x[k] + g_y_0_xzzz_xyyz[k];

                g_y_0_xxzzz_yzz[k] = -g_y_0_xzzz_yzz[k] * cd_x[k] + g_y_0_xzzz_xyzz[k];

                g_y_0_xxzzz_zzz[k] = -g_y_0_xzzz_zzz[k] * cd_x[k] + g_y_0_xzzz_xzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 100);

            auto g_y_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 101);

            auto g_y_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 102);

            auto g_y_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 103);

            auto g_y_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 104);

            auto g_y_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 105);

            auto g_y_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 106);

            auto g_y_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 107);

            auto g_y_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 108);

            auto g_y_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 109);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyy_xxx, g_y_0_xyyyy_xxy, g_y_0_xyyyy_xxz, g_y_0_xyyyy_xyy, g_y_0_xyyyy_xyz, g_y_0_xyyyy_xzz, g_y_0_xyyyy_yyy, g_y_0_xyyyy_yyz, g_y_0_xyyyy_yzz, g_y_0_xyyyy_zzz, g_y_0_yyyy_xxx, g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxy, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xyy, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_yyy, g_y_0_yyyy_yyz, g_y_0_yyyy_yzz, g_y_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_xxx[k] = -g_y_0_yyyy_xxx[k] * cd_x[k] + g_y_0_yyyy_xxxx[k];

                g_y_0_xyyyy_xxy[k] = -g_y_0_yyyy_xxy[k] * cd_x[k] + g_y_0_yyyy_xxxy[k];

                g_y_0_xyyyy_xxz[k] = -g_y_0_yyyy_xxz[k] * cd_x[k] + g_y_0_yyyy_xxxz[k];

                g_y_0_xyyyy_xyy[k] = -g_y_0_yyyy_xyy[k] * cd_x[k] + g_y_0_yyyy_xxyy[k];

                g_y_0_xyyyy_xyz[k] = -g_y_0_yyyy_xyz[k] * cd_x[k] + g_y_0_yyyy_xxyz[k];

                g_y_0_xyyyy_xzz[k] = -g_y_0_yyyy_xzz[k] * cd_x[k] + g_y_0_yyyy_xxzz[k];

                g_y_0_xyyyy_yyy[k] = -g_y_0_yyyy_yyy[k] * cd_x[k] + g_y_0_yyyy_xyyy[k];

                g_y_0_xyyyy_yyz[k] = -g_y_0_yyyy_yyz[k] * cd_x[k] + g_y_0_yyyy_xyyz[k];

                g_y_0_xyyyy_yzz[k] = -g_y_0_yyyy_yzz[k] * cd_x[k] + g_y_0_yyyy_xyzz[k];

                g_y_0_xyyyy_zzz[k] = -g_y_0_yyyy_zzz[k] * cd_x[k] + g_y_0_yyyy_xzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 110);

            auto g_y_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 111);

            auto g_y_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 112);

            auto g_y_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 113);

            auto g_y_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 114);

            auto g_y_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 115);

            auto g_y_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 116);

            auto g_y_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 117);

            auto g_y_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 118);

            auto g_y_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyz_xxx, g_y_0_xyyyz_xxy, g_y_0_xyyyz_xxz, g_y_0_xyyyz_xyy, g_y_0_xyyyz_xyz, g_y_0_xyyyz_xzz, g_y_0_xyyyz_yyy, g_y_0_xyyyz_yyz, g_y_0_xyyyz_yzz, g_y_0_xyyyz_zzz, g_y_0_yyyz_xxx, g_y_0_yyyz_xxxx, g_y_0_yyyz_xxxy, g_y_0_yyyz_xxxz, g_y_0_yyyz_xxy, g_y_0_yyyz_xxyy, g_y_0_yyyz_xxyz, g_y_0_yyyz_xxz, g_y_0_yyyz_xxzz, g_y_0_yyyz_xyy, g_y_0_yyyz_xyyy, g_y_0_yyyz_xyyz, g_y_0_yyyz_xyz, g_y_0_yyyz_xyzz, g_y_0_yyyz_xzz, g_y_0_yyyz_xzzz, g_y_0_yyyz_yyy, g_y_0_yyyz_yyz, g_y_0_yyyz_yzz, g_y_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_xxx[k] = -g_y_0_yyyz_xxx[k] * cd_x[k] + g_y_0_yyyz_xxxx[k];

                g_y_0_xyyyz_xxy[k] = -g_y_0_yyyz_xxy[k] * cd_x[k] + g_y_0_yyyz_xxxy[k];

                g_y_0_xyyyz_xxz[k] = -g_y_0_yyyz_xxz[k] * cd_x[k] + g_y_0_yyyz_xxxz[k];

                g_y_0_xyyyz_xyy[k] = -g_y_0_yyyz_xyy[k] * cd_x[k] + g_y_0_yyyz_xxyy[k];

                g_y_0_xyyyz_xyz[k] = -g_y_0_yyyz_xyz[k] * cd_x[k] + g_y_0_yyyz_xxyz[k];

                g_y_0_xyyyz_xzz[k] = -g_y_0_yyyz_xzz[k] * cd_x[k] + g_y_0_yyyz_xxzz[k];

                g_y_0_xyyyz_yyy[k] = -g_y_0_yyyz_yyy[k] * cd_x[k] + g_y_0_yyyz_xyyy[k];

                g_y_0_xyyyz_yyz[k] = -g_y_0_yyyz_yyz[k] * cd_x[k] + g_y_0_yyyz_xyyz[k];

                g_y_0_xyyyz_yzz[k] = -g_y_0_yyyz_yzz[k] * cd_x[k] + g_y_0_yyyz_xyzz[k];

                g_y_0_xyyyz_zzz[k] = -g_y_0_yyyz_zzz[k] * cd_x[k] + g_y_0_yyyz_xzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 120);

            auto g_y_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 121);

            auto g_y_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 122);

            auto g_y_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 123);

            auto g_y_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 124);

            auto g_y_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 125);

            auto g_y_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 126);

            auto g_y_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 127);

            auto g_y_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 128);

            auto g_y_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 129);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzz_xxx, g_y_0_xyyzz_xxy, g_y_0_xyyzz_xxz, g_y_0_xyyzz_xyy, g_y_0_xyyzz_xyz, g_y_0_xyyzz_xzz, g_y_0_xyyzz_yyy, g_y_0_xyyzz_yyz, g_y_0_xyyzz_yzz, g_y_0_xyyzz_zzz, g_y_0_yyzz_xxx, g_y_0_yyzz_xxxx, g_y_0_yyzz_xxxy, g_y_0_yyzz_xxxz, g_y_0_yyzz_xxy, g_y_0_yyzz_xxyy, g_y_0_yyzz_xxyz, g_y_0_yyzz_xxz, g_y_0_yyzz_xxzz, g_y_0_yyzz_xyy, g_y_0_yyzz_xyyy, g_y_0_yyzz_xyyz, g_y_0_yyzz_xyz, g_y_0_yyzz_xyzz, g_y_0_yyzz_xzz, g_y_0_yyzz_xzzz, g_y_0_yyzz_yyy, g_y_0_yyzz_yyz, g_y_0_yyzz_yzz, g_y_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_xxx[k] = -g_y_0_yyzz_xxx[k] * cd_x[k] + g_y_0_yyzz_xxxx[k];

                g_y_0_xyyzz_xxy[k] = -g_y_0_yyzz_xxy[k] * cd_x[k] + g_y_0_yyzz_xxxy[k];

                g_y_0_xyyzz_xxz[k] = -g_y_0_yyzz_xxz[k] * cd_x[k] + g_y_0_yyzz_xxxz[k];

                g_y_0_xyyzz_xyy[k] = -g_y_0_yyzz_xyy[k] * cd_x[k] + g_y_0_yyzz_xxyy[k];

                g_y_0_xyyzz_xyz[k] = -g_y_0_yyzz_xyz[k] * cd_x[k] + g_y_0_yyzz_xxyz[k];

                g_y_0_xyyzz_xzz[k] = -g_y_0_yyzz_xzz[k] * cd_x[k] + g_y_0_yyzz_xxzz[k];

                g_y_0_xyyzz_yyy[k] = -g_y_0_yyzz_yyy[k] * cd_x[k] + g_y_0_yyzz_xyyy[k];

                g_y_0_xyyzz_yyz[k] = -g_y_0_yyzz_yyz[k] * cd_x[k] + g_y_0_yyzz_xyyz[k];

                g_y_0_xyyzz_yzz[k] = -g_y_0_yyzz_yzz[k] * cd_x[k] + g_y_0_yyzz_xyzz[k];

                g_y_0_xyyzz_zzz[k] = -g_y_0_yyzz_zzz[k] * cd_x[k] + g_y_0_yyzz_xzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 130);

            auto g_y_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 131);

            auto g_y_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 132);

            auto g_y_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 133);

            auto g_y_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 134);

            auto g_y_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 135);

            auto g_y_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 136);

            auto g_y_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 137);

            auto g_y_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 138);

            auto g_y_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzz_xxx, g_y_0_xyzzz_xxy, g_y_0_xyzzz_xxz, g_y_0_xyzzz_xyy, g_y_0_xyzzz_xyz, g_y_0_xyzzz_xzz, g_y_0_xyzzz_yyy, g_y_0_xyzzz_yyz, g_y_0_xyzzz_yzz, g_y_0_xyzzz_zzz, g_y_0_yzzz_xxx, g_y_0_yzzz_xxxx, g_y_0_yzzz_xxxy, g_y_0_yzzz_xxxz, g_y_0_yzzz_xxy, g_y_0_yzzz_xxyy, g_y_0_yzzz_xxyz, g_y_0_yzzz_xxz, g_y_0_yzzz_xxzz, g_y_0_yzzz_xyy, g_y_0_yzzz_xyyy, g_y_0_yzzz_xyyz, g_y_0_yzzz_xyz, g_y_0_yzzz_xyzz, g_y_0_yzzz_xzz, g_y_0_yzzz_xzzz, g_y_0_yzzz_yyy, g_y_0_yzzz_yyz, g_y_0_yzzz_yzz, g_y_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_xxx[k] = -g_y_0_yzzz_xxx[k] * cd_x[k] + g_y_0_yzzz_xxxx[k];

                g_y_0_xyzzz_xxy[k] = -g_y_0_yzzz_xxy[k] * cd_x[k] + g_y_0_yzzz_xxxy[k];

                g_y_0_xyzzz_xxz[k] = -g_y_0_yzzz_xxz[k] * cd_x[k] + g_y_0_yzzz_xxxz[k];

                g_y_0_xyzzz_xyy[k] = -g_y_0_yzzz_xyy[k] * cd_x[k] + g_y_0_yzzz_xxyy[k];

                g_y_0_xyzzz_xyz[k] = -g_y_0_yzzz_xyz[k] * cd_x[k] + g_y_0_yzzz_xxyz[k];

                g_y_0_xyzzz_xzz[k] = -g_y_0_yzzz_xzz[k] * cd_x[k] + g_y_0_yzzz_xxzz[k];

                g_y_0_xyzzz_yyy[k] = -g_y_0_yzzz_yyy[k] * cd_x[k] + g_y_0_yzzz_xyyy[k];

                g_y_0_xyzzz_yyz[k] = -g_y_0_yzzz_yyz[k] * cd_x[k] + g_y_0_yzzz_xyyz[k];

                g_y_0_xyzzz_yzz[k] = -g_y_0_yzzz_yzz[k] * cd_x[k] + g_y_0_yzzz_xyzz[k];

                g_y_0_xyzzz_zzz[k] = -g_y_0_yzzz_zzz[k] * cd_x[k] + g_y_0_yzzz_xzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 140);

            auto g_y_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 141);

            auto g_y_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 142);

            auto g_y_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 143);

            auto g_y_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 144);

            auto g_y_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 145);

            auto g_y_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 146);

            auto g_y_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 147);

            auto g_y_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 148);

            auto g_y_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzz_xxx, g_y_0_xzzzz_xxy, g_y_0_xzzzz_xxz, g_y_0_xzzzz_xyy, g_y_0_xzzzz_xyz, g_y_0_xzzzz_xzz, g_y_0_xzzzz_yyy, g_y_0_xzzzz_yyz, g_y_0_xzzzz_yzz, g_y_0_xzzzz_zzz, g_y_0_zzzz_xxx, g_y_0_zzzz_xxxx, g_y_0_zzzz_xxxy, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxy, g_y_0_zzzz_xxyy, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xyy, g_y_0_zzzz_xyyy, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_yyy, g_y_0_zzzz_yyz, g_y_0_zzzz_yzz, g_y_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_xxx[k] = -g_y_0_zzzz_xxx[k] * cd_x[k] + g_y_0_zzzz_xxxx[k];

                g_y_0_xzzzz_xxy[k] = -g_y_0_zzzz_xxy[k] * cd_x[k] + g_y_0_zzzz_xxxy[k];

                g_y_0_xzzzz_xxz[k] = -g_y_0_zzzz_xxz[k] * cd_x[k] + g_y_0_zzzz_xxxz[k];

                g_y_0_xzzzz_xyy[k] = -g_y_0_zzzz_xyy[k] * cd_x[k] + g_y_0_zzzz_xxyy[k];

                g_y_0_xzzzz_xyz[k] = -g_y_0_zzzz_xyz[k] * cd_x[k] + g_y_0_zzzz_xxyz[k];

                g_y_0_xzzzz_xzz[k] = -g_y_0_zzzz_xzz[k] * cd_x[k] + g_y_0_zzzz_xxzz[k];

                g_y_0_xzzzz_yyy[k] = -g_y_0_zzzz_yyy[k] * cd_x[k] + g_y_0_zzzz_xyyy[k];

                g_y_0_xzzzz_yyz[k] = -g_y_0_zzzz_yyz[k] * cd_x[k] + g_y_0_zzzz_xyyz[k];

                g_y_0_xzzzz_yzz[k] = -g_y_0_zzzz_yzz[k] * cd_x[k] + g_y_0_zzzz_xyzz[k];

                g_y_0_xzzzz_zzz[k] = -g_y_0_zzzz_zzz[k] * cd_x[k] + g_y_0_zzzz_xzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 150);

            auto g_y_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 151);

            auto g_y_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 152);

            auto g_y_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 153);

            auto g_y_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 154);

            auto g_y_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 155);

            auto g_y_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 156);

            auto g_y_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 157);

            auto g_y_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 158);

            auto g_y_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 159);

            #pragma omp simd aligned(cd_y, g_y_0_yyyy_xxx, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxy, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxz, g_y_0_yyyy_xyy, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xzz, g_y_0_yyyy_yyy, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_zzz, g_y_0_yyyyy_xxx, g_y_0_yyyyy_xxy, g_y_0_yyyyy_xxz, g_y_0_yyyyy_xyy, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xzz, g_y_0_yyyyy_yyy, g_y_0_yyyyy_yyz, g_y_0_yyyyy_yzz, g_y_0_yyyyy_zzz, g_yyyy_xxx, g_yyyy_xxy, g_yyyy_xxz, g_yyyy_xyy, g_yyyy_xyz, g_yyyy_xzz, g_yyyy_yyy, g_yyyy_yyz, g_yyyy_yzz, g_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_xxx[k] = -g_yyyy_xxx[k] - g_y_0_yyyy_xxx[k] * cd_y[k] + g_y_0_yyyy_xxxy[k];

                g_y_0_yyyyy_xxy[k] = -g_yyyy_xxy[k] - g_y_0_yyyy_xxy[k] * cd_y[k] + g_y_0_yyyy_xxyy[k];

                g_y_0_yyyyy_xxz[k] = -g_yyyy_xxz[k] - g_y_0_yyyy_xxz[k] * cd_y[k] + g_y_0_yyyy_xxyz[k];

                g_y_0_yyyyy_xyy[k] = -g_yyyy_xyy[k] - g_y_0_yyyy_xyy[k] * cd_y[k] + g_y_0_yyyy_xyyy[k];

                g_y_0_yyyyy_xyz[k] = -g_yyyy_xyz[k] - g_y_0_yyyy_xyz[k] * cd_y[k] + g_y_0_yyyy_xyyz[k];

                g_y_0_yyyyy_xzz[k] = -g_yyyy_xzz[k] - g_y_0_yyyy_xzz[k] * cd_y[k] + g_y_0_yyyy_xyzz[k];

                g_y_0_yyyyy_yyy[k] = -g_yyyy_yyy[k] - g_y_0_yyyy_yyy[k] * cd_y[k] + g_y_0_yyyy_yyyy[k];

                g_y_0_yyyyy_yyz[k] = -g_yyyy_yyz[k] - g_y_0_yyyy_yyz[k] * cd_y[k] + g_y_0_yyyy_yyyz[k];

                g_y_0_yyyyy_yzz[k] = -g_yyyy_yzz[k] - g_y_0_yyyy_yzz[k] * cd_y[k] + g_y_0_yyyy_yyzz[k];

                g_y_0_yyyyy_zzz[k] = -g_yyyy_zzz[k] - g_y_0_yyyy_zzz[k] * cd_y[k] + g_y_0_yyyy_yzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 160);

            auto g_y_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 161);

            auto g_y_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 162);

            auto g_y_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 163);

            auto g_y_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 164);

            auto g_y_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 165);

            auto g_y_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 166);

            auto g_y_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 167);

            auto g_y_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 168);

            auto g_y_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 169);

            #pragma omp simd aligned(cd_z, g_y_0_yyyy_xxx, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxy, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xyy, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_yyy, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_zzz, g_y_0_yyyy_zzzz, g_y_0_yyyyz_xxx, g_y_0_yyyyz_xxy, g_y_0_yyyyz_xxz, g_y_0_yyyyz_xyy, g_y_0_yyyyz_xyz, g_y_0_yyyyz_xzz, g_y_0_yyyyz_yyy, g_y_0_yyyyz_yyz, g_y_0_yyyyz_yzz, g_y_0_yyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_xxx[k] = -g_y_0_yyyy_xxx[k] * cd_z[k] + g_y_0_yyyy_xxxz[k];

                g_y_0_yyyyz_xxy[k] = -g_y_0_yyyy_xxy[k] * cd_z[k] + g_y_0_yyyy_xxyz[k];

                g_y_0_yyyyz_xxz[k] = -g_y_0_yyyy_xxz[k] * cd_z[k] + g_y_0_yyyy_xxzz[k];

                g_y_0_yyyyz_xyy[k] = -g_y_0_yyyy_xyy[k] * cd_z[k] + g_y_0_yyyy_xyyz[k];

                g_y_0_yyyyz_xyz[k] = -g_y_0_yyyy_xyz[k] * cd_z[k] + g_y_0_yyyy_xyzz[k];

                g_y_0_yyyyz_xzz[k] = -g_y_0_yyyy_xzz[k] * cd_z[k] + g_y_0_yyyy_xzzz[k];

                g_y_0_yyyyz_yyy[k] = -g_y_0_yyyy_yyy[k] * cd_z[k] + g_y_0_yyyy_yyyz[k];

                g_y_0_yyyyz_yyz[k] = -g_y_0_yyyy_yyz[k] * cd_z[k] + g_y_0_yyyy_yyzz[k];

                g_y_0_yyyyz_yzz[k] = -g_y_0_yyyy_yzz[k] * cd_z[k] + g_y_0_yyyy_yzzz[k];

                g_y_0_yyyyz_zzz[k] = -g_y_0_yyyy_zzz[k] * cd_z[k] + g_y_0_yyyy_zzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 170);

            auto g_y_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 171);

            auto g_y_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 172);

            auto g_y_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 173);

            auto g_y_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 174);

            auto g_y_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 175);

            auto g_y_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 176);

            auto g_y_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 177);

            auto g_y_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 178);

            auto g_y_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_z, g_y_0_yyyz_xxx, g_y_0_yyyz_xxxz, g_y_0_yyyz_xxy, g_y_0_yyyz_xxyz, g_y_0_yyyz_xxz, g_y_0_yyyz_xxzz, g_y_0_yyyz_xyy, g_y_0_yyyz_xyyz, g_y_0_yyyz_xyz, g_y_0_yyyz_xyzz, g_y_0_yyyz_xzz, g_y_0_yyyz_xzzz, g_y_0_yyyz_yyy, g_y_0_yyyz_yyyz, g_y_0_yyyz_yyz, g_y_0_yyyz_yyzz, g_y_0_yyyz_yzz, g_y_0_yyyz_yzzz, g_y_0_yyyz_zzz, g_y_0_yyyz_zzzz, g_y_0_yyyzz_xxx, g_y_0_yyyzz_xxy, g_y_0_yyyzz_xxz, g_y_0_yyyzz_xyy, g_y_0_yyyzz_xyz, g_y_0_yyyzz_xzz, g_y_0_yyyzz_yyy, g_y_0_yyyzz_yyz, g_y_0_yyyzz_yzz, g_y_0_yyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_xxx[k] = -g_y_0_yyyz_xxx[k] * cd_z[k] + g_y_0_yyyz_xxxz[k];

                g_y_0_yyyzz_xxy[k] = -g_y_0_yyyz_xxy[k] * cd_z[k] + g_y_0_yyyz_xxyz[k];

                g_y_0_yyyzz_xxz[k] = -g_y_0_yyyz_xxz[k] * cd_z[k] + g_y_0_yyyz_xxzz[k];

                g_y_0_yyyzz_xyy[k] = -g_y_0_yyyz_xyy[k] * cd_z[k] + g_y_0_yyyz_xyyz[k];

                g_y_0_yyyzz_xyz[k] = -g_y_0_yyyz_xyz[k] * cd_z[k] + g_y_0_yyyz_xyzz[k];

                g_y_0_yyyzz_xzz[k] = -g_y_0_yyyz_xzz[k] * cd_z[k] + g_y_0_yyyz_xzzz[k];

                g_y_0_yyyzz_yyy[k] = -g_y_0_yyyz_yyy[k] * cd_z[k] + g_y_0_yyyz_yyyz[k];

                g_y_0_yyyzz_yyz[k] = -g_y_0_yyyz_yyz[k] * cd_z[k] + g_y_0_yyyz_yyzz[k];

                g_y_0_yyyzz_yzz[k] = -g_y_0_yyyz_yzz[k] * cd_z[k] + g_y_0_yyyz_yzzz[k];

                g_y_0_yyyzz_zzz[k] = -g_y_0_yyyz_zzz[k] * cd_z[k] + g_y_0_yyyz_zzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 180);

            auto g_y_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 181);

            auto g_y_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 182);

            auto g_y_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 183);

            auto g_y_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 184);

            auto g_y_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 185);

            auto g_y_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 186);

            auto g_y_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 187);

            auto g_y_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 188);

            auto g_y_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 189);

            #pragma omp simd aligned(cd_z, g_y_0_yyzz_xxx, g_y_0_yyzz_xxxz, g_y_0_yyzz_xxy, g_y_0_yyzz_xxyz, g_y_0_yyzz_xxz, g_y_0_yyzz_xxzz, g_y_0_yyzz_xyy, g_y_0_yyzz_xyyz, g_y_0_yyzz_xyz, g_y_0_yyzz_xyzz, g_y_0_yyzz_xzz, g_y_0_yyzz_xzzz, g_y_0_yyzz_yyy, g_y_0_yyzz_yyyz, g_y_0_yyzz_yyz, g_y_0_yyzz_yyzz, g_y_0_yyzz_yzz, g_y_0_yyzz_yzzz, g_y_0_yyzz_zzz, g_y_0_yyzz_zzzz, g_y_0_yyzzz_xxx, g_y_0_yyzzz_xxy, g_y_0_yyzzz_xxz, g_y_0_yyzzz_xyy, g_y_0_yyzzz_xyz, g_y_0_yyzzz_xzz, g_y_0_yyzzz_yyy, g_y_0_yyzzz_yyz, g_y_0_yyzzz_yzz, g_y_0_yyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_xxx[k] = -g_y_0_yyzz_xxx[k] * cd_z[k] + g_y_0_yyzz_xxxz[k];

                g_y_0_yyzzz_xxy[k] = -g_y_0_yyzz_xxy[k] * cd_z[k] + g_y_0_yyzz_xxyz[k];

                g_y_0_yyzzz_xxz[k] = -g_y_0_yyzz_xxz[k] * cd_z[k] + g_y_0_yyzz_xxzz[k];

                g_y_0_yyzzz_xyy[k] = -g_y_0_yyzz_xyy[k] * cd_z[k] + g_y_0_yyzz_xyyz[k];

                g_y_0_yyzzz_xyz[k] = -g_y_0_yyzz_xyz[k] * cd_z[k] + g_y_0_yyzz_xyzz[k];

                g_y_0_yyzzz_xzz[k] = -g_y_0_yyzz_xzz[k] * cd_z[k] + g_y_0_yyzz_xzzz[k];

                g_y_0_yyzzz_yyy[k] = -g_y_0_yyzz_yyy[k] * cd_z[k] + g_y_0_yyzz_yyyz[k];

                g_y_0_yyzzz_yyz[k] = -g_y_0_yyzz_yyz[k] * cd_z[k] + g_y_0_yyzz_yyzz[k];

                g_y_0_yyzzz_yzz[k] = -g_y_0_yyzz_yzz[k] * cd_z[k] + g_y_0_yyzz_yzzz[k];

                g_y_0_yyzzz_zzz[k] = -g_y_0_yyzz_zzz[k] * cd_z[k] + g_y_0_yyzz_zzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 190);

            auto g_y_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 191);

            auto g_y_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 192);

            auto g_y_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 193);

            auto g_y_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 194);

            auto g_y_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 195);

            auto g_y_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 196);

            auto g_y_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 197);

            auto g_y_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 198);

            auto g_y_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 199);

            #pragma omp simd aligned(cd_z, g_y_0_yzzz_xxx, g_y_0_yzzz_xxxz, g_y_0_yzzz_xxy, g_y_0_yzzz_xxyz, g_y_0_yzzz_xxz, g_y_0_yzzz_xxzz, g_y_0_yzzz_xyy, g_y_0_yzzz_xyyz, g_y_0_yzzz_xyz, g_y_0_yzzz_xyzz, g_y_0_yzzz_xzz, g_y_0_yzzz_xzzz, g_y_0_yzzz_yyy, g_y_0_yzzz_yyyz, g_y_0_yzzz_yyz, g_y_0_yzzz_yyzz, g_y_0_yzzz_yzz, g_y_0_yzzz_yzzz, g_y_0_yzzz_zzz, g_y_0_yzzz_zzzz, g_y_0_yzzzz_xxx, g_y_0_yzzzz_xxy, g_y_0_yzzzz_xxz, g_y_0_yzzzz_xyy, g_y_0_yzzzz_xyz, g_y_0_yzzzz_xzz, g_y_0_yzzzz_yyy, g_y_0_yzzzz_yyz, g_y_0_yzzzz_yzz, g_y_0_yzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_xxx[k] = -g_y_0_yzzz_xxx[k] * cd_z[k] + g_y_0_yzzz_xxxz[k];

                g_y_0_yzzzz_xxy[k] = -g_y_0_yzzz_xxy[k] * cd_z[k] + g_y_0_yzzz_xxyz[k];

                g_y_0_yzzzz_xxz[k] = -g_y_0_yzzz_xxz[k] * cd_z[k] + g_y_0_yzzz_xxzz[k];

                g_y_0_yzzzz_xyy[k] = -g_y_0_yzzz_xyy[k] * cd_z[k] + g_y_0_yzzz_xyyz[k];

                g_y_0_yzzzz_xyz[k] = -g_y_0_yzzz_xyz[k] * cd_z[k] + g_y_0_yzzz_xyzz[k];

                g_y_0_yzzzz_xzz[k] = -g_y_0_yzzz_xzz[k] * cd_z[k] + g_y_0_yzzz_xzzz[k];

                g_y_0_yzzzz_yyy[k] = -g_y_0_yzzz_yyy[k] * cd_z[k] + g_y_0_yzzz_yyyz[k];

                g_y_0_yzzzz_yyz[k] = -g_y_0_yzzz_yyz[k] * cd_z[k] + g_y_0_yzzz_yyzz[k];

                g_y_0_yzzzz_yzz[k] = -g_y_0_yzzz_yzz[k] * cd_z[k] + g_y_0_yzzz_yzzz[k];

                g_y_0_yzzzz_zzz[k] = -g_y_0_yzzz_zzz[k] * cd_z[k] + g_y_0_yzzz_zzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 200);

            auto g_y_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 201);

            auto g_y_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 202);

            auto g_y_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 203);

            auto g_y_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 204);

            auto g_y_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 205);

            auto g_y_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 206);

            auto g_y_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 207);

            auto g_y_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 208);

            auto g_y_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_z, g_y_0_zzzz_xxx, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxy, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xyy, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_yyy, g_y_0_zzzz_yyyz, g_y_0_zzzz_yyz, g_y_0_zzzz_yyzz, g_y_0_zzzz_yzz, g_y_0_zzzz_yzzz, g_y_0_zzzz_zzz, g_y_0_zzzz_zzzz, g_y_0_zzzzz_xxx, g_y_0_zzzzz_xxy, g_y_0_zzzzz_xxz, g_y_0_zzzzz_xyy, g_y_0_zzzzz_xyz, g_y_0_zzzzz_xzz, g_y_0_zzzzz_yyy, g_y_0_zzzzz_yyz, g_y_0_zzzzz_yzz, g_y_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_xxx[k] = -g_y_0_zzzz_xxx[k] * cd_z[k] + g_y_0_zzzz_xxxz[k];

                g_y_0_zzzzz_xxy[k] = -g_y_0_zzzz_xxy[k] * cd_z[k] + g_y_0_zzzz_xxyz[k];

                g_y_0_zzzzz_xxz[k] = -g_y_0_zzzz_xxz[k] * cd_z[k] + g_y_0_zzzz_xxzz[k];

                g_y_0_zzzzz_xyy[k] = -g_y_0_zzzz_xyy[k] * cd_z[k] + g_y_0_zzzz_xyyz[k];

                g_y_0_zzzzz_xyz[k] = -g_y_0_zzzz_xyz[k] * cd_z[k] + g_y_0_zzzz_xyzz[k];

                g_y_0_zzzzz_xzz[k] = -g_y_0_zzzz_xzz[k] * cd_z[k] + g_y_0_zzzz_xzzz[k];

                g_y_0_zzzzz_yyy[k] = -g_y_0_zzzz_yyy[k] * cd_z[k] + g_y_0_zzzz_yyyz[k];

                g_y_0_zzzzz_yyz[k] = -g_y_0_zzzz_yyz[k] * cd_z[k] + g_y_0_zzzz_yyzz[k];

                g_y_0_zzzzz_yzz[k] = -g_y_0_zzzz_yzz[k] * cd_z[k] + g_y_0_zzzz_yzzz[k];

                g_y_0_zzzzz_zzz[k] = -g_y_0_zzzz_zzz[k] * cd_z[k] + g_y_0_zzzz_zzzz[k];
            }
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 5);

            auto g_z_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 6);

            auto g_z_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 7);

            auto g_z_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 8);

            auto g_z_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_z_0_xxxx_xxx, g_z_0_xxxx_xxxx, g_z_0_xxxx_xxxy, g_z_0_xxxx_xxxz, g_z_0_xxxx_xxy, g_z_0_xxxx_xxyy, g_z_0_xxxx_xxyz, g_z_0_xxxx_xxz, g_z_0_xxxx_xxzz, g_z_0_xxxx_xyy, g_z_0_xxxx_xyyy, g_z_0_xxxx_xyyz, g_z_0_xxxx_xyz, g_z_0_xxxx_xyzz, g_z_0_xxxx_xzz, g_z_0_xxxx_xzzz, g_z_0_xxxx_yyy, g_z_0_xxxx_yyz, g_z_0_xxxx_yzz, g_z_0_xxxx_zzz, g_z_0_xxxxx_xxx, g_z_0_xxxxx_xxy, g_z_0_xxxxx_xxz, g_z_0_xxxxx_xyy, g_z_0_xxxxx_xyz, g_z_0_xxxxx_xzz, g_z_0_xxxxx_yyy, g_z_0_xxxxx_yyz, g_z_0_xxxxx_yzz, g_z_0_xxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_xxx[k] = -g_z_0_xxxx_xxx[k] * cd_x[k] + g_z_0_xxxx_xxxx[k];

                g_z_0_xxxxx_xxy[k] = -g_z_0_xxxx_xxy[k] * cd_x[k] + g_z_0_xxxx_xxxy[k];

                g_z_0_xxxxx_xxz[k] = -g_z_0_xxxx_xxz[k] * cd_x[k] + g_z_0_xxxx_xxxz[k];

                g_z_0_xxxxx_xyy[k] = -g_z_0_xxxx_xyy[k] * cd_x[k] + g_z_0_xxxx_xxyy[k];

                g_z_0_xxxxx_xyz[k] = -g_z_0_xxxx_xyz[k] * cd_x[k] + g_z_0_xxxx_xxyz[k];

                g_z_0_xxxxx_xzz[k] = -g_z_0_xxxx_xzz[k] * cd_x[k] + g_z_0_xxxx_xxzz[k];

                g_z_0_xxxxx_yyy[k] = -g_z_0_xxxx_yyy[k] * cd_x[k] + g_z_0_xxxx_xyyy[k];

                g_z_0_xxxxx_yyz[k] = -g_z_0_xxxx_yyz[k] * cd_x[k] + g_z_0_xxxx_xyyz[k];

                g_z_0_xxxxx_yzz[k] = -g_z_0_xxxx_yzz[k] * cd_x[k] + g_z_0_xxxx_xyzz[k];

                g_z_0_xxxxx_zzz[k] = -g_z_0_xxxx_zzz[k] * cd_x[k] + g_z_0_xxxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 10);

            auto g_z_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 11);

            auto g_z_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 12);

            auto g_z_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 13);

            auto g_z_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 14);

            auto g_z_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 15);

            auto g_z_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 16);

            auto g_z_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 17);

            auto g_z_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 18);

            auto g_z_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxy_xxx, g_z_0_xxxxy_xxy, g_z_0_xxxxy_xxz, g_z_0_xxxxy_xyy, g_z_0_xxxxy_xyz, g_z_0_xxxxy_xzz, g_z_0_xxxxy_yyy, g_z_0_xxxxy_yyz, g_z_0_xxxxy_yzz, g_z_0_xxxxy_zzz, g_z_0_xxxy_xxx, g_z_0_xxxy_xxxx, g_z_0_xxxy_xxxy, g_z_0_xxxy_xxxz, g_z_0_xxxy_xxy, g_z_0_xxxy_xxyy, g_z_0_xxxy_xxyz, g_z_0_xxxy_xxz, g_z_0_xxxy_xxzz, g_z_0_xxxy_xyy, g_z_0_xxxy_xyyy, g_z_0_xxxy_xyyz, g_z_0_xxxy_xyz, g_z_0_xxxy_xyzz, g_z_0_xxxy_xzz, g_z_0_xxxy_xzzz, g_z_0_xxxy_yyy, g_z_0_xxxy_yyz, g_z_0_xxxy_yzz, g_z_0_xxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_xxx[k] = -g_z_0_xxxy_xxx[k] * cd_x[k] + g_z_0_xxxy_xxxx[k];

                g_z_0_xxxxy_xxy[k] = -g_z_0_xxxy_xxy[k] * cd_x[k] + g_z_0_xxxy_xxxy[k];

                g_z_0_xxxxy_xxz[k] = -g_z_0_xxxy_xxz[k] * cd_x[k] + g_z_0_xxxy_xxxz[k];

                g_z_0_xxxxy_xyy[k] = -g_z_0_xxxy_xyy[k] * cd_x[k] + g_z_0_xxxy_xxyy[k];

                g_z_0_xxxxy_xyz[k] = -g_z_0_xxxy_xyz[k] * cd_x[k] + g_z_0_xxxy_xxyz[k];

                g_z_0_xxxxy_xzz[k] = -g_z_0_xxxy_xzz[k] * cd_x[k] + g_z_0_xxxy_xxzz[k];

                g_z_0_xxxxy_yyy[k] = -g_z_0_xxxy_yyy[k] * cd_x[k] + g_z_0_xxxy_xyyy[k];

                g_z_0_xxxxy_yyz[k] = -g_z_0_xxxy_yyz[k] * cd_x[k] + g_z_0_xxxy_xyyz[k];

                g_z_0_xxxxy_yzz[k] = -g_z_0_xxxy_yzz[k] * cd_x[k] + g_z_0_xxxy_xyzz[k];

                g_z_0_xxxxy_zzz[k] = -g_z_0_xxxy_zzz[k] * cd_x[k] + g_z_0_xxxy_xzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 20);

            auto g_z_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 21);

            auto g_z_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 22);

            auto g_z_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 23);

            auto g_z_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 24);

            auto g_z_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 25);

            auto g_z_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 26);

            auto g_z_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 27);

            auto g_z_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 28);

            auto g_z_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxz_xxx, g_z_0_xxxxz_xxy, g_z_0_xxxxz_xxz, g_z_0_xxxxz_xyy, g_z_0_xxxxz_xyz, g_z_0_xxxxz_xzz, g_z_0_xxxxz_yyy, g_z_0_xxxxz_yyz, g_z_0_xxxxz_yzz, g_z_0_xxxxz_zzz, g_z_0_xxxz_xxx, g_z_0_xxxz_xxxx, g_z_0_xxxz_xxxy, g_z_0_xxxz_xxxz, g_z_0_xxxz_xxy, g_z_0_xxxz_xxyy, g_z_0_xxxz_xxyz, g_z_0_xxxz_xxz, g_z_0_xxxz_xxzz, g_z_0_xxxz_xyy, g_z_0_xxxz_xyyy, g_z_0_xxxz_xyyz, g_z_0_xxxz_xyz, g_z_0_xxxz_xyzz, g_z_0_xxxz_xzz, g_z_0_xxxz_xzzz, g_z_0_xxxz_yyy, g_z_0_xxxz_yyz, g_z_0_xxxz_yzz, g_z_0_xxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_xxx[k] = -g_z_0_xxxz_xxx[k] * cd_x[k] + g_z_0_xxxz_xxxx[k];

                g_z_0_xxxxz_xxy[k] = -g_z_0_xxxz_xxy[k] * cd_x[k] + g_z_0_xxxz_xxxy[k];

                g_z_0_xxxxz_xxz[k] = -g_z_0_xxxz_xxz[k] * cd_x[k] + g_z_0_xxxz_xxxz[k];

                g_z_0_xxxxz_xyy[k] = -g_z_0_xxxz_xyy[k] * cd_x[k] + g_z_0_xxxz_xxyy[k];

                g_z_0_xxxxz_xyz[k] = -g_z_0_xxxz_xyz[k] * cd_x[k] + g_z_0_xxxz_xxyz[k];

                g_z_0_xxxxz_xzz[k] = -g_z_0_xxxz_xzz[k] * cd_x[k] + g_z_0_xxxz_xxzz[k];

                g_z_0_xxxxz_yyy[k] = -g_z_0_xxxz_yyy[k] * cd_x[k] + g_z_0_xxxz_xyyy[k];

                g_z_0_xxxxz_yyz[k] = -g_z_0_xxxz_yyz[k] * cd_x[k] + g_z_0_xxxz_xyyz[k];

                g_z_0_xxxxz_yzz[k] = -g_z_0_xxxz_yzz[k] * cd_x[k] + g_z_0_xxxz_xyzz[k];

                g_z_0_xxxxz_zzz[k] = -g_z_0_xxxz_zzz[k] * cd_x[k] + g_z_0_xxxz_xzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 30);

            auto g_z_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 31);

            auto g_z_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 32);

            auto g_z_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 33);

            auto g_z_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 34);

            auto g_z_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 35);

            auto g_z_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 36);

            auto g_z_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 37);

            auto g_z_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 38);

            auto g_z_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 39);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyy_xxx, g_z_0_xxxyy_xxy, g_z_0_xxxyy_xxz, g_z_0_xxxyy_xyy, g_z_0_xxxyy_xyz, g_z_0_xxxyy_xzz, g_z_0_xxxyy_yyy, g_z_0_xxxyy_yyz, g_z_0_xxxyy_yzz, g_z_0_xxxyy_zzz, g_z_0_xxyy_xxx, g_z_0_xxyy_xxxx, g_z_0_xxyy_xxxy, g_z_0_xxyy_xxxz, g_z_0_xxyy_xxy, g_z_0_xxyy_xxyy, g_z_0_xxyy_xxyz, g_z_0_xxyy_xxz, g_z_0_xxyy_xxzz, g_z_0_xxyy_xyy, g_z_0_xxyy_xyyy, g_z_0_xxyy_xyyz, g_z_0_xxyy_xyz, g_z_0_xxyy_xyzz, g_z_0_xxyy_xzz, g_z_0_xxyy_xzzz, g_z_0_xxyy_yyy, g_z_0_xxyy_yyz, g_z_0_xxyy_yzz, g_z_0_xxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_xxx[k] = -g_z_0_xxyy_xxx[k] * cd_x[k] + g_z_0_xxyy_xxxx[k];

                g_z_0_xxxyy_xxy[k] = -g_z_0_xxyy_xxy[k] * cd_x[k] + g_z_0_xxyy_xxxy[k];

                g_z_0_xxxyy_xxz[k] = -g_z_0_xxyy_xxz[k] * cd_x[k] + g_z_0_xxyy_xxxz[k];

                g_z_0_xxxyy_xyy[k] = -g_z_0_xxyy_xyy[k] * cd_x[k] + g_z_0_xxyy_xxyy[k];

                g_z_0_xxxyy_xyz[k] = -g_z_0_xxyy_xyz[k] * cd_x[k] + g_z_0_xxyy_xxyz[k];

                g_z_0_xxxyy_xzz[k] = -g_z_0_xxyy_xzz[k] * cd_x[k] + g_z_0_xxyy_xxzz[k];

                g_z_0_xxxyy_yyy[k] = -g_z_0_xxyy_yyy[k] * cd_x[k] + g_z_0_xxyy_xyyy[k];

                g_z_0_xxxyy_yyz[k] = -g_z_0_xxyy_yyz[k] * cd_x[k] + g_z_0_xxyy_xyyz[k];

                g_z_0_xxxyy_yzz[k] = -g_z_0_xxyy_yzz[k] * cd_x[k] + g_z_0_xxyy_xyzz[k];

                g_z_0_xxxyy_zzz[k] = -g_z_0_xxyy_zzz[k] * cd_x[k] + g_z_0_xxyy_xzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 40);

            auto g_z_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 41);

            auto g_z_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 42);

            auto g_z_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 43);

            auto g_z_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 44);

            auto g_z_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 45);

            auto g_z_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 46);

            auto g_z_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 47);

            auto g_z_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 48);

            auto g_z_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 49);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyz_xxx, g_z_0_xxxyz_xxy, g_z_0_xxxyz_xxz, g_z_0_xxxyz_xyy, g_z_0_xxxyz_xyz, g_z_0_xxxyz_xzz, g_z_0_xxxyz_yyy, g_z_0_xxxyz_yyz, g_z_0_xxxyz_yzz, g_z_0_xxxyz_zzz, g_z_0_xxyz_xxx, g_z_0_xxyz_xxxx, g_z_0_xxyz_xxxy, g_z_0_xxyz_xxxz, g_z_0_xxyz_xxy, g_z_0_xxyz_xxyy, g_z_0_xxyz_xxyz, g_z_0_xxyz_xxz, g_z_0_xxyz_xxzz, g_z_0_xxyz_xyy, g_z_0_xxyz_xyyy, g_z_0_xxyz_xyyz, g_z_0_xxyz_xyz, g_z_0_xxyz_xyzz, g_z_0_xxyz_xzz, g_z_0_xxyz_xzzz, g_z_0_xxyz_yyy, g_z_0_xxyz_yyz, g_z_0_xxyz_yzz, g_z_0_xxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_xxx[k] = -g_z_0_xxyz_xxx[k] * cd_x[k] + g_z_0_xxyz_xxxx[k];

                g_z_0_xxxyz_xxy[k] = -g_z_0_xxyz_xxy[k] * cd_x[k] + g_z_0_xxyz_xxxy[k];

                g_z_0_xxxyz_xxz[k] = -g_z_0_xxyz_xxz[k] * cd_x[k] + g_z_0_xxyz_xxxz[k];

                g_z_0_xxxyz_xyy[k] = -g_z_0_xxyz_xyy[k] * cd_x[k] + g_z_0_xxyz_xxyy[k];

                g_z_0_xxxyz_xyz[k] = -g_z_0_xxyz_xyz[k] * cd_x[k] + g_z_0_xxyz_xxyz[k];

                g_z_0_xxxyz_xzz[k] = -g_z_0_xxyz_xzz[k] * cd_x[k] + g_z_0_xxyz_xxzz[k];

                g_z_0_xxxyz_yyy[k] = -g_z_0_xxyz_yyy[k] * cd_x[k] + g_z_0_xxyz_xyyy[k];

                g_z_0_xxxyz_yyz[k] = -g_z_0_xxyz_yyz[k] * cd_x[k] + g_z_0_xxyz_xyyz[k];

                g_z_0_xxxyz_yzz[k] = -g_z_0_xxyz_yzz[k] * cd_x[k] + g_z_0_xxyz_xyzz[k];

                g_z_0_xxxyz_zzz[k] = -g_z_0_xxyz_zzz[k] * cd_x[k] + g_z_0_xxyz_xzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 50);

            auto g_z_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 51);

            auto g_z_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 52);

            auto g_z_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 53);

            auto g_z_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 54);

            auto g_z_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 55);

            auto g_z_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 56);

            auto g_z_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 57);

            auto g_z_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 58);

            auto g_z_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzz_xxx, g_z_0_xxxzz_xxy, g_z_0_xxxzz_xxz, g_z_0_xxxzz_xyy, g_z_0_xxxzz_xyz, g_z_0_xxxzz_xzz, g_z_0_xxxzz_yyy, g_z_0_xxxzz_yyz, g_z_0_xxxzz_yzz, g_z_0_xxxzz_zzz, g_z_0_xxzz_xxx, g_z_0_xxzz_xxxx, g_z_0_xxzz_xxxy, g_z_0_xxzz_xxxz, g_z_0_xxzz_xxy, g_z_0_xxzz_xxyy, g_z_0_xxzz_xxyz, g_z_0_xxzz_xxz, g_z_0_xxzz_xxzz, g_z_0_xxzz_xyy, g_z_0_xxzz_xyyy, g_z_0_xxzz_xyyz, g_z_0_xxzz_xyz, g_z_0_xxzz_xyzz, g_z_0_xxzz_xzz, g_z_0_xxzz_xzzz, g_z_0_xxzz_yyy, g_z_0_xxzz_yyz, g_z_0_xxzz_yzz, g_z_0_xxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_xxx[k] = -g_z_0_xxzz_xxx[k] * cd_x[k] + g_z_0_xxzz_xxxx[k];

                g_z_0_xxxzz_xxy[k] = -g_z_0_xxzz_xxy[k] * cd_x[k] + g_z_0_xxzz_xxxy[k];

                g_z_0_xxxzz_xxz[k] = -g_z_0_xxzz_xxz[k] * cd_x[k] + g_z_0_xxzz_xxxz[k];

                g_z_0_xxxzz_xyy[k] = -g_z_0_xxzz_xyy[k] * cd_x[k] + g_z_0_xxzz_xxyy[k];

                g_z_0_xxxzz_xyz[k] = -g_z_0_xxzz_xyz[k] * cd_x[k] + g_z_0_xxzz_xxyz[k];

                g_z_0_xxxzz_xzz[k] = -g_z_0_xxzz_xzz[k] * cd_x[k] + g_z_0_xxzz_xxzz[k];

                g_z_0_xxxzz_yyy[k] = -g_z_0_xxzz_yyy[k] * cd_x[k] + g_z_0_xxzz_xyyy[k];

                g_z_0_xxxzz_yyz[k] = -g_z_0_xxzz_yyz[k] * cd_x[k] + g_z_0_xxzz_xyyz[k];

                g_z_0_xxxzz_yzz[k] = -g_z_0_xxzz_yzz[k] * cd_x[k] + g_z_0_xxzz_xyzz[k];

                g_z_0_xxxzz_zzz[k] = -g_z_0_xxzz_zzz[k] * cd_x[k] + g_z_0_xxzz_xzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 60);

            auto g_z_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 61);

            auto g_z_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 62);

            auto g_z_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 63);

            auto g_z_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 64);

            auto g_z_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 65);

            auto g_z_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 66);

            auto g_z_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 67);

            auto g_z_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 68);

            auto g_z_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 69);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyy_xxx, g_z_0_xxyyy_xxy, g_z_0_xxyyy_xxz, g_z_0_xxyyy_xyy, g_z_0_xxyyy_xyz, g_z_0_xxyyy_xzz, g_z_0_xxyyy_yyy, g_z_0_xxyyy_yyz, g_z_0_xxyyy_yzz, g_z_0_xxyyy_zzz, g_z_0_xyyy_xxx, g_z_0_xyyy_xxxx, g_z_0_xyyy_xxxy, g_z_0_xyyy_xxxz, g_z_0_xyyy_xxy, g_z_0_xyyy_xxyy, g_z_0_xyyy_xxyz, g_z_0_xyyy_xxz, g_z_0_xyyy_xxzz, g_z_0_xyyy_xyy, g_z_0_xyyy_xyyy, g_z_0_xyyy_xyyz, g_z_0_xyyy_xyz, g_z_0_xyyy_xyzz, g_z_0_xyyy_xzz, g_z_0_xyyy_xzzz, g_z_0_xyyy_yyy, g_z_0_xyyy_yyz, g_z_0_xyyy_yzz, g_z_0_xyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_xxx[k] = -g_z_0_xyyy_xxx[k] * cd_x[k] + g_z_0_xyyy_xxxx[k];

                g_z_0_xxyyy_xxy[k] = -g_z_0_xyyy_xxy[k] * cd_x[k] + g_z_0_xyyy_xxxy[k];

                g_z_0_xxyyy_xxz[k] = -g_z_0_xyyy_xxz[k] * cd_x[k] + g_z_0_xyyy_xxxz[k];

                g_z_0_xxyyy_xyy[k] = -g_z_0_xyyy_xyy[k] * cd_x[k] + g_z_0_xyyy_xxyy[k];

                g_z_0_xxyyy_xyz[k] = -g_z_0_xyyy_xyz[k] * cd_x[k] + g_z_0_xyyy_xxyz[k];

                g_z_0_xxyyy_xzz[k] = -g_z_0_xyyy_xzz[k] * cd_x[k] + g_z_0_xyyy_xxzz[k];

                g_z_0_xxyyy_yyy[k] = -g_z_0_xyyy_yyy[k] * cd_x[k] + g_z_0_xyyy_xyyy[k];

                g_z_0_xxyyy_yyz[k] = -g_z_0_xyyy_yyz[k] * cd_x[k] + g_z_0_xyyy_xyyz[k];

                g_z_0_xxyyy_yzz[k] = -g_z_0_xyyy_yzz[k] * cd_x[k] + g_z_0_xyyy_xyzz[k];

                g_z_0_xxyyy_zzz[k] = -g_z_0_xyyy_zzz[k] * cd_x[k] + g_z_0_xyyy_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 70);

            auto g_z_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 71);

            auto g_z_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 72);

            auto g_z_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 73);

            auto g_z_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 74);

            auto g_z_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 75);

            auto g_z_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 76);

            auto g_z_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 77);

            auto g_z_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 78);

            auto g_z_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 79);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyz_xxx, g_z_0_xxyyz_xxy, g_z_0_xxyyz_xxz, g_z_0_xxyyz_xyy, g_z_0_xxyyz_xyz, g_z_0_xxyyz_xzz, g_z_0_xxyyz_yyy, g_z_0_xxyyz_yyz, g_z_0_xxyyz_yzz, g_z_0_xxyyz_zzz, g_z_0_xyyz_xxx, g_z_0_xyyz_xxxx, g_z_0_xyyz_xxxy, g_z_0_xyyz_xxxz, g_z_0_xyyz_xxy, g_z_0_xyyz_xxyy, g_z_0_xyyz_xxyz, g_z_0_xyyz_xxz, g_z_0_xyyz_xxzz, g_z_0_xyyz_xyy, g_z_0_xyyz_xyyy, g_z_0_xyyz_xyyz, g_z_0_xyyz_xyz, g_z_0_xyyz_xyzz, g_z_0_xyyz_xzz, g_z_0_xyyz_xzzz, g_z_0_xyyz_yyy, g_z_0_xyyz_yyz, g_z_0_xyyz_yzz, g_z_0_xyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_xxx[k] = -g_z_0_xyyz_xxx[k] * cd_x[k] + g_z_0_xyyz_xxxx[k];

                g_z_0_xxyyz_xxy[k] = -g_z_0_xyyz_xxy[k] * cd_x[k] + g_z_0_xyyz_xxxy[k];

                g_z_0_xxyyz_xxz[k] = -g_z_0_xyyz_xxz[k] * cd_x[k] + g_z_0_xyyz_xxxz[k];

                g_z_0_xxyyz_xyy[k] = -g_z_0_xyyz_xyy[k] * cd_x[k] + g_z_0_xyyz_xxyy[k];

                g_z_0_xxyyz_xyz[k] = -g_z_0_xyyz_xyz[k] * cd_x[k] + g_z_0_xyyz_xxyz[k];

                g_z_0_xxyyz_xzz[k] = -g_z_0_xyyz_xzz[k] * cd_x[k] + g_z_0_xyyz_xxzz[k];

                g_z_0_xxyyz_yyy[k] = -g_z_0_xyyz_yyy[k] * cd_x[k] + g_z_0_xyyz_xyyy[k];

                g_z_0_xxyyz_yyz[k] = -g_z_0_xyyz_yyz[k] * cd_x[k] + g_z_0_xyyz_xyyz[k];

                g_z_0_xxyyz_yzz[k] = -g_z_0_xyyz_yzz[k] * cd_x[k] + g_z_0_xyyz_xyzz[k];

                g_z_0_xxyyz_zzz[k] = -g_z_0_xyyz_zzz[k] * cd_x[k] + g_z_0_xyyz_xzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 80);

            auto g_z_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 81);

            auto g_z_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 82);

            auto g_z_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 83);

            auto g_z_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 84);

            auto g_z_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 85);

            auto g_z_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 86);

            auto g_z_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 87);

            auto g_z_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 88);

            auto g_z_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzz_xxx, g_z_0_xxyzz_xxy, g_z_0_xxyzz_xxz, g_z_0_xxyzz_xyy, g_z_0_xxyzz_xyz, g_z_0_xxyzz_xzz, g_z_0_xxyzz_yyy, g_z_0_xxyzz_yyz, g_z_0_xxyzz_yzz, g_z_0_xxyzz_zzz, g_z_0_xyzz_xxx, g_z_0_xyzz_xxxx, g_z_0_xyzz_xxxy, g_z_0_xyzz_xxxz, g_z_0_xyzz_xxy, g_z_0_xyzz_xxyy, g_z_0_xyzz_xxyz, g_z_0_xyzz_xxz, g_z_0_xyzz_xxzz, g_z_0_xyzz_xyy, g_z_0_xyzz_xyyy, g_z_0_xyzz_xyyz, g_z_0_xyzz_xyz, g_z_0_xyzz_xyzz, g_z_0_xyzz_xzz, g_z_0_xyzz_xzzz, g_z_0_xyzz_yyy, g_z_0_xyzz_yyz, g_z_0_xyzz_yzz, g_z_0_xyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_xxx[k] = -g_z_0_xyzz_xxx[k] * cd_x[k] + g_z_0_xyzz_xxxx[k];

                g_z_0_xxyzz_xxy[k] = -g_z_0_xyzz_xxy[k] * cd_x[k] + g_z_0_xyzz_xxxy[k];

                g_z_0_xxyzz_xxz[k] = -g_z_0_xyzz_xxz[k] * cd_x[k] + g_z_0_xyzz_xxxz[k];

                g_z_0_xxyzz_xyy[k] = -g_z_0_xyzz_xyy[k] * cd_x[k] + g_z_0_xyzz_xxyy[k];

                g_z_0_xxyzz_xyz[k] = -g_z_0_xyzz_xyz[k] * cd_x[k] + g_z_0_xyzz_xxyz[k];

                g_z_0_xxyzz_xzz[k] = -g_z_0_xyzz_xzz[k] * cd_x[k] + g_z_0_xyzz_xxzz[k];

                g_z_0_xxyzz_yyy[k] = -g_z_0_xyzz_yyy[k] * cd_x[k] + g_z_0_xyzz_xyyy[k];

                g_z_0_xxyzz_yyz[k] = -g_z_0_xyzz_yyz[k] * cd_x[k] + g_z_0_xyzz_xyyz[k];

                g_z_0_xxyzz_yzz[k] = -g_z_0_xyzz_yzz[k] * cd_x[k] + g_z_0_xyzz_xyzz[k];

                g_z_0_xxyzz_zzz[k] = -g_z_0_xyzz_zzz[k] * cd_x[k] + g_z_0_xyzz_xzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 90);

            auto g_z_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 91);

            auto g_z_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 92);

            auto g_z_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 93);

            auto g_z_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 94);

            auto g_z_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 95);

            auto g_z_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 96);

            auto g_z_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 97);

            auto g_z_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 98);

            auto g_z_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 99);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzz_xxx, g_z_0_xxzzz_xxy, g_z_0_xxzzz_xxz, g_z_0_xxzzz_xyy, g_z_0_xxzzz_xyz, g_z_0_xxzzz_xzz, g_z_0_xxzzz_yyy, g_z_0_xxzzz_yyz, g_z_0_xxzzz_yzz, g_z_0_xxzzz_zzz, g_z_0_xzzz_xxx, g_z_0_xzzz_xxxx, g_z_0_xzzz_xxxy, g_z_0_xzzz_xxxz, g_z_0_xzzz_xxy, g_z_0_xzzz_xxyy, g_z_0_xzzz_xxyz, g_z_0_xzzz_xxz, g_z_0_xzzz_xxzz, g_z_0_xzzz_xyy, g_z_0_xzzz_xyyy, g_z_0_xzzz_xyyz, g_z_0_xzzz_xyz, g_z_0_xzzz_xyzz, g_z_0_xzzz_xzz, g_z_0_xzzz_xzzz, g_z_0_xzzz_yyy, g_z_0_xzzz_yyz, g_z_0_xzzz_yzz, g_z_0_xzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_xxx[k] = -g_z_0_xzzz_xxx[k] * cd_x[k] + g_z_0_xzzz_xxxx[k];

                g_z_0_xxzzz_xxy[k] = -g_z_0_xzzz_xxy[k] * cd_x[k] + g_z_0_xzzz_xxxy[k];

                g_z_0_xxzzz_xxz[k] = -g_z_0_xzzz_xxz[k] * cd_x[k] + g_z_0_xzzz_xxxz[k];

                g_z_0_xxzzz_xyy[k] = -g_z_0_xzzz_xyy[k] * cd_x[k] + g_z_0_xzzz_xxyy[k];

                g_z_0_xxzzz_xyz[k] = -g_z_0_xzzz_xyz[k] * cd_x[k] + g_z_0_xzzz_xxyz[k];

                g_z_0_xxzzz_xzz[k] = -g_z_0_xzzz_xzz[k] * cd_x[k] + g_z_0_xzzz_xxzz[k];

                g_z_0_xxzzz_yyy[k] = -g_z_0_xzzz_yyy[k] * cd_x[k] + g_z_0_xzzz_xyyy[k];

                g_z_0_xxzzz_yyz[k] = -g_z_0_xzzz_yyz[k] * cd_x[k] + g_z_0_xzzz_xyyz[k];

                g_z_0_xxzzz_yzz[k] = -g_z_0_xzzz_yzz[k] * cd_x[k] + g_z_0_xzzz_xyzz[k];

                g_z_0_xxzzz_zzz[k] = -g_z_0_xzzz_zzz[k] * cd_x[k] + g_z_0_xzzz_xzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 100);

            auto g_z_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 101);

            auto g_z_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 102);

            auto g_z_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 103);

            auto g_z_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 104);

            auto g_z_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 105);

            auto g_z_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 106);

            auto g_z_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 107);

            auto g_z_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 108);

            auto g_z_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 109);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyy_xxx, g_z_0_xyyyy_xxy, g_z_0_xyyyy_xxz, g_z_0_xyyyy_xyy, g_z_0_xyyyy_xyz, g_z_0_xyyyy_xzz, g_z_0_xyyyy_yyy, g_z_0_xyyyy_yyz, g_z_0_xyyyy_yzz, g_z_0_xyyyy_zzz, g_z_0_yyyy_xxx, g_z_0_yyyy_xxxx, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxxz, g_z_0_yyyy_xxy, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxz, g_z_0_yyyy_xxzz, g_z_0_yyyy_xyy, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xzz, g_z_0_yyyy_xzzz, g_z_0_yyyy_yyy, g_z_0_yyyy_yyz, g_z_0_yyyy_yzz, g_z_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_xxx[k] = -g_z_0_yyyy_xxx[k] * cd_x[k] + g_z_0_yyyy_xxxx[k];

                g_z_0_xyyyy_xxy[k] = -g_z_0_yyyy_xxy[k] * cd_x[k] + g_z_0_yyyy_xxxy[k];

                g_z_0_xyyyy_xxz[k] = -g_z_0_yyyy_xxz[k] * cd_x[k] + g_z_0_yyyy_xxxz[k];

                g_z_0_xyyyy_xyy[k] = -g_z_0_yyyy_xyy[k] * cd_x[k] + g_z_0_yyyy_xxyy[k];

                g_z_0_xyyyy_xyz[k] = -g_z_0_yyyy_xyz[k] * cd_x[k] + g_z_0_yyyy_xxyz[k];

                g_z_0_xyyyy_xzz[k] = -g_z_0_yyyy_xzz[k] * cd_x[k] + g_z_0_yyyy_xxzz[k];

                g_z_0_xyyyy_yyy[k] = -g_z_0_yyyy_yyy[k] * cd_x[k] + g_z_0_yyyy_xyyy[k];

                g_z_0_xyyyy_yyz[k] = -g_z_0_yyyy_yyz[k] * cd_x[k] + g_z_0_yyyy_xyyz[k];

                g_z_0_xyyyy_yzz[k] = -g_z_0_yyyy_yzz[k] * cd_x[k] + g_z_0_yyyy_xyzz[k];

                g_z_0_xyyyy_zzz[k] = -g_z_0_yyyy_zzz[k] * cd_x[k] + g_z_0_yyyy_xzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 110);

            auto g_z_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 111);

            auto g_z_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 112);

            auto g_z_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 113);

            auto g_z_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 114);

            auto g_z_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 115);

            auto g_z_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 116);

            auto g_z_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 117);

            auto g_z_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 118);

            auto g_z_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyz_xxx, g_z_0_xyyyz_xxy, g_z_0_xyyyz_xxz, g_z_0_xyyyz_xyy, g_z_0_xyyyz_xyz, g_z_0_xyyyz_xzz, g_z_0_xyyyz_yyy, g_z_0_xyyyz_yyz, g_z_0_xyyyz_yzz, g_z_0_xyyyz_zzz, g_z_0_yyyz_xxx, g_z_0_yyyz_xxxx, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxxz, g_z_0_yyyz_xxy, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxz, g_z_0_yyyz_xxzz, g_z_0_yyyz_xyy, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xzz, g_z_0_yyyz_xzzz, g_z_0_yyyz_yyy, g_z_0_yyyz_yyz, g_z_0_yyyz_yzz, g_z_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_xxx[k] = -g_z_0_yyyz_xxx[k] * cd_x[k] + g_z_0_yyyz_xxxx[k];

                g_z_0_xyyyz_xxy[k] = -g_z_0_yyyz_xxy[k] * cd_x[k] + g_z_0_yyyz_xxxy[k];

                g_z_0_xyyyz_xxz[k] = -g_z_0_yyyz_xxz[k] * cd_x[k] + g_z_0_yyyz_xxxz[k];

                g_z_0_xyyyz_xyy[k] = -g_z_0_yyyz_xyy[k] * cd_x[k] + g_z_0_yyyz_xxyy[k];

                g_z_0_xyyyz_xyz[k] = -g_z_0_yyyz_xyz[k] * cd_x[k] + g_z_0_yyyz_xxyz[k];

                g_z_0_xyyyz_xzz[k] = -g_z_0_yyyz_xzz[k] * cd_x[k] + g_z_0_yyyz_xxzz[k];

                g_z_0_xyyyz_yyy[k] = -g_z_0_yyyz_yyy[k] * cd_x[k] + g_z_0_yyyz_xyyy[k];

                g_z_0_xyyyz_yyz[k] = -g_z_0_yyyz_yyz[k] * cd_x[k] + g_z_0_yyyz_xyyz[k];

                g_z_0_xyyyz_yzz[k] = -g_z_0_yyyz_yzz[k] * cd_x[k] + g_z_0_yyyz_xyzz[k];

                g_z_0_xyyyz_zzz[k] = -g_z_0_yyyz_zzz[k] * cd_x[k] + g_z_0_yyyz_xzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 120);

            auto g_z_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 121);

            auto g_z_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 122);

            auto g_z_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 123);

            auto g_z_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 124);

            auto g_z_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 125);

            auto g_z_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 126);

            auto g_z_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 127);

            auto g_z_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 128);

            auto g_z_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 129);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzz_xxx, g_z_0_xyyzz_xxy, g_z_0_xyyzz_xxz, g_z_0_xyyzz_xyy, g_z_0_xyyzz_xyz, g_z_0_xyyzz_xzz, g_z_0_xyyzz_yyy, g_z_0_xyyzz_yyz, g_z_0_xyyzz_yzz, g_z_0_xyyzz_zzz, g_z_0_yyzz_xxx, g_z_0_yyzz_xxxx, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxxz, g_z_0_yyzz_xxy, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxz, g_z_0_yyzz_xxzz, g_z_0_yyzz_xyy, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xzz, g_z_0_yyzz_xzzz, g_z_0_yyzz_yyy, g_z_0_yyzz_yyz, g_z_0_yyzz_yzz, g_z_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_xxx[k] = -g_z_0_yyzz_xxx[k] * cd_x[k] + g_z_0_yyzz_xxxx[k];

                g_z_0_xyyzz_xxy[k] = -g_z_0_yyzz_xxy[k] * cd_x[k] + g_z_0_yyzz_xxxy[k];

                g_z_0_xyyzz_xxz[k] = -g_z_0_yyzz_xxz[k] * cd_x[k] + g_z_0_yyzz_xxxz[k];

                g_z_0_xyyzz_xyy[k] = -g_z_0_yyzz_xyy[k] * cd_x[k] + g_z_0_yyzz_xxyy[k];

                g_z_0_xyyzz_xyz[k] = -g_z_0_yyzz_xyz[k] * cd_x[k] + g_z_0_yyzz_xxyz[k];

                g_z_0_xyyzz_xzz[k] = -g_z_0_yyzz_xzz[k] * cd_x[k] + g_z_0_yyzz_xxzz[k];

                g_z_0_xyyzz_yyy[k] = -g_z_0_yyzz_yyy[k] * cd_x[k] + g_z_0_yyzz_xyyy[k];

                g_z_0_xyyzz_yyz[k] = -g_z_0_yyzz_yyz[k] * cd_x[k] + g_z_0_yyzz_xyyz[k];

                g_z_0_xyyzz_yzz[k] = -g_z_0_yyzz_yzz[k] * cd_x[k] + g_z_0_yyzz_xyzz[k];

                g_z_0_xyyzz_zzz[k] = -g_z_0_yyzz_zzz[k] * cd_x[k] + g_z_0_yyzz_xzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 130);

            auto g_z_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 131);

            auto g_z_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 132);

            auto g_z_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 133);

            auto g_z_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 134);

            auto g_z_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 135);

            auto g_z_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 136);

            auto g_z_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 137);

            auto g_z_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 138);

            auto g_z_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzz_xxx, g_z_0_xyzzz_xxy, g_z_0_xyzzz_xxz, g_z_0_xyzzz_xyy, g_z_0_xyzzz_xyz, g_z_0_xyzzz_xzz, g_z_0_xyzzz_yyy, g_z_0_xyzzz_yyz, g_z_0_xyzzz_yzz, g_z_0_xyzzz_zzz, g_z_0_yzzz_xxx, g_z_0_yzzz_xxxx, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxxz, g_z_0_yzzz_xxy, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxz, g_z_0_yzzz_xxzz, g_z_0_yzzz_xyy, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xzz, g_z_0_yzzz_xzzz, g_z_0_yzzz_yyy, g_z_0_yzzz_yyz, g_z_0_yzzz_yzz, g_z_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_xxx[k] = -g_z_0_yzzz_xxx[k] * cd_x[k] + g_z_0_yzzz_xxxx[k];

                g_z_0_xyzzz_xxy[k] = -g_z_0_yzzz_xxy[k] * cd_x[k] + g_z_0_yzzz_xxxy[k];

                g_z_0_xyzzz_xxz[k] = -g_z_0_yzzz_xxz[k] * cd_x[k] + g_z_0_yzzz_xxxz[k];

                g_z_0_xyzzz_xyy[k] = -g_z_0_yzzz_xyy[k] * cd_x[k] + g_z_0_yzzz_xxyy[k];

                g_z_0_xyzzz_xyz[k] = -g_z_0_yzzz_xyz[k] * cd_x[k] + g_z_0_yzzz_xxyz[k];

                g_z_0_xyzzz_xzz[k] = -g_z_0_yzzz_xzz[k] * cd_x[k] + g_z_0_yzzz_xxzz[k];

                g_z_0_xyzzz_yyy[k] = -g_z_0_yzzz_yyy[k] * cd_x[k] + g_z_0_yzzz_xyyy[k];

                g_z_0_xyzzz_yyz[k] = -g_z_0_yzzz_yyz[k] * cd_x[k] + g_z_0_yzzz_xyyz[k];

                g_z_0_xyzzz_yzz[k] = -g_z_0_yzzz_yzz[k] * cd_x[k] + g_z_0_yzzz_xyzz[k];

                g_z_0_xyzzz_zzz[k] = -g_z_0_yzzz_zzz[k] * cd_x[k] + g_z_0_yzzz_xzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 140);

            auto g_z_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 141);

            auto g_z_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 142);

            auto g_z_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 143);

            auto g_z_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 144);

            auto g_z_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 145);

            auto g_z_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 146);

            auto g_z_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 147);

            auto g_z_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 148);

            auto g_z_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzz_xxx, g_z_0_xzzzz_xxy, g_z_0_xzzzz_xxz, g_z_0_xzzzz_xyy, g_z_0_xzzzz_xyz, g_z_0_xzzzz_xzz, g_z_0_xzzzz_yyy, g_z_0_xzzzz_yyz, g_z_0_xzzzz_yzz, g_z_0_xzzzz_zzz, g_z_0_zzzz_xxx, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxy, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyy, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyy, g_z_0_zzzz_yyz, g_z_0_zzzz_yzz, g_z_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_xxx[k] = -g_z_0_zzzz_xxx[k] * cd_x[k] + g_z_0_zzzz_xxxx[k];

                g_z_0_xzzzz_xxy[k] = -g_z_0_zzzz_xxy[k] * cd_x[k] + g_z_0_zzzz_xxxy[k];

                g_z_0_xzzzz_xxz[k] = -g_z_0_zzzz_xxz[k] * cd_x[k] + g_z_0_zzzz_xxxz[k];

                g_z_0_xzzzz_xyy[k] = -g_z_0_zzzz_xyy[k] * cd_x[k] + g_z_0_zzzz_xxyy[k];

                g_z_0_xzzzz_xyz[k] = -g_z_0_zzzz_xyz[k] * cd_x[k] + g_z_0_zzzz_xxyz[k];

                g_z_0_xzzzz_xzz[k] = -g_z_0_zzzz_xzz[k] * cd_x[k] + g_z_0_zzzz_xxzz[k];

                g_z_0_xzzzz_yyy[k] = -g_z_0_zzzz_yyy[k] * cd_x[k] + g_z_0_zzzz_xyyy[k];

                g_z_0_xzzzz_yyz[k] = -g_z_0_zzzz_yyz[k] * cd_x[k] + g_z_0_zzzz_xyyz[k];

                g_z_0_xzzzz_yzz[k] = -g_z_0_zzzz_yzz[k] * cd_x[k] + g_z_0_zzzz_xyzz[k];

                g_z_0_xzzzz_zzz[k] = -g_z_0_zzzz_zzz[k] * cd_x[k] + g_z_0_zzzz_xzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 150);

            auto g_z_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 151);

            auto g_z_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 152);

            auto g_z_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 153);

            auto g_z_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 154);

            auto g_z_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 155);

            auto g_z_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 156);

            auto g_z_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 157);

            auto g_z_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 158);

            auto g_z_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 159);

            #pragma omp simd aligned(cd_y, g_z_0_yyyy_xxx, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxy, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxz, g_z_0_yyyy_xyy, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xzz, g_z_0_yyyy_yyy, g_z_0_yyyy_yyyy, g_z_0_yyyy_yyyz, g_z_0_yyyy_yyz, g_z_0_yyyy_yyzz, g_z_0_yyyy_yzz, g_z_0_yyyy_yzzz, g_z_0_yyyy_zzz, g_z_0_yyyyy_xxx, g_z_0_yyyyy_xxy, g_z_0_yyyyy_xxz, g_z_0_yyyyy_xyy, g_z_0_yyyyy_xyz, g_z_0_yyyyy_xzz, g_z_0_yyyyy_yyy, g_z_0_yyyyy_yyz, g_z_0_yyyyy_yzz, g_z_0_yyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_xxx[k] = -g_z_0_yyyy_xxx[k] * cd_y[k] + g_z_0_yyyy_xxxy[k];

                g_z_0_yyyyy_xxy[k] = -g_z_0_yyyy_xxy[k] * cd_y[k] + g_z_0_yyyy_xxyy[k];

                g_z_0_yyyyy_xxz[k] = -g_z_0_yyyy_xxz[k] * cd_y[k] + g_z_0_yyyy_xxyz[k];

                g_z_0_yyyyy_xyy[k] = -g_z_0_yyyy_xyy[k] * cd_y[k] + g_z_0_yyyy_xyyy[k];

                g_z_0_yyyyy_xyz[k] = -g_z_0_yyyy_xyz[k] * cd_y[k] + g_z_0_yyyy_xyyz[k];

                g_z_0_yyyyy_xzz[k] = -g_z_0_yyyy_xzz[k] * cd_y[k] + g_z_0_yyyy_xyzz[k];

                g_z_0_yyyyy_yyy[k] = -g_z_0_yyyy_yyy[k] * cd_y[k] + g_z_0_yyyy_yyyy[k];

                g_z_0_yyyyy_yyz[k] = -g_z_0_yyyy_yyz[k] * cd_y[k] + g_z_0_yyyy_yyyz[k];

                g_z_0_yyyyy_yzz[k] = -g_z_0_yyyy_yzz[k] * cd_y[k] + g_z_0_yyyy_yyzz[k];

                g_z_0_yyyyy_zzz[k] = -g_z_0_yyyy_zzz[k] * cd_y[k] + g_z_0_yyyy_yzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 160);

            auto g_z_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 161);

            auto g_z_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 162);

            auto g_z_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 163);

            auto g_z_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 164);

            auto g_z_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 165);

            auto g_z_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 166);

            auto g_z_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 167);

            auto g_z_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 168);

            auto g_z_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 169);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyz_xxx, g_z_0_yyyyz_xxy, g_z_0_yyyyz_xxz, g_z_0_yyyyz_xyy, g_z_0_yyyyz_xyz, g_z_0_yyyyz_xzz, g_z_0_yyyyz_yyy, g_z_0_yyyyz_yyz, g_z_0_yyyyz_yzz, g_z_0_yyyyz_zzz, g_z_0_yyyz_xxx, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxy, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxz, g_z_0_yyyz_xyy, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xzz, g_z_0_yyyz_yyy, g_z_0_yyyz_yyyy, g_z_0_yyyz_yyyz, g_z_0_yyyz_yyz, g_z_0_yyyz_yyzz, g_z_0_yyyz_yzz, g_z_0_yyyz_yzzz, g_z_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_xxx[k] = -g_z_0_yyyz_xxx[k] * cd_y[k] + g_z_0_yyyz_xxxy[k];

                g_z_0_yyyyz_xxy[k] = -g_z_0_yyyz_xxy[k] * cd_y[k] + g_z_0_yyyz_xxyy[k];

                g_z_0_yyyyz_xxz[k] = -g_z_0_yyyz_xxz[k] * cd_y[k] + g_z_0_yyyz_xxyz[k];

                g_z_0_yyyyz_xyy[k] = -g_z_0_yyyz_xyy[k] * cd_y[k] + g_z_0_yyyz_xyyy[k];

                g_z_0_yyyyz_xyz[k] = -g_z_0_yyyz_xyz[k] * cd_y[k] + g_z_0_yyyz_xyyz[k];

                g_z_0_yyyyz_xzz[k] = -g_z_0_yyyz_xzz[k] * cd_y[k] + g_z_0_yyyz_xyzz[k];

                g_z_0_yyyyz_yyy[k] = -g_z_0_yyyz_yyy[k] * cd_y[k] + g_z_0_yyyz_yyyy[k];

                g_z_0_yyyyz_yyz[k] = -g_z_0_yyyz_yyz[k] * cd_y[k] + g_z_0_yyyz_yyyz[k];

                g_z_0_yyyyz_yzz[k] = -g_z_0_yyyz_yzz[k] * cd_y[k] + g_z_0_yyyz_yyzz[k];

                g_z_0_yyyyz_zzz[k] = -g_z_0_yyyz_zzz[k] * cd_y[k] + g_z_0_yyyz_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 170);

            auto g_z_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 171);

            auto g_z_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 172);

            auto g_z_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 173);

            auto g_z_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 174);

            auto g_z_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 175);

            auto g_z_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 176);

            auto g_z_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 177);

            auto g_z_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 178);

            auto g_z_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzz_xxx, g_z_0_yyyzz_xxy, g_z_0_yyyzz_xxz, g_z_0_yyyzz_xyy, g_z_0_yyyzz_xyz, g_z_0_yyyzz_xzz, g_z_0_yyyzz_yyy, g_z_0_yyyzz_yyz, g_z_0_yyyzz_yzz, g_z_0_yyyzz_zzz, g_z_0_yyzz_xxx, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxy, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxz, g_z_0_yyzz_xyy, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xzz, g_z_0_yyzz_yyy, g_z_0_yyzz_yyyy, g_z_0_yyzz_yyyz, g_z_0_yyzz_yyz, g_z_0_yyzz_yyzz, g_z_0_yyzz_yzz, g_z_0_yyzz_yzzz, g_z_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_xxx[k] = -g_z_0_yyzz_xxx[k] * cd_y[k] + g_z_0_yyzz_xxxy[k];

                g_z_0_yyyzz_xxy[k] = -g_z_0_yyzz_xxy[k] * cd_y[k] + g_z_0_yyzz_xxyy[k];

                g_z_0_yyyzz_xxz[k] = -g_z_0_yyzz_xxz[k] * cd_y[k] + g_z_0_yyzz_xxyz[k];

                g_z_0_yyyzz_xyy[k] = -g_z_0_yyzz_xyy[k] * cd_y[k] + g_z_0_yyzz_xyyy[k];

                g_z_0_yyyzz_xyz[k] = -g_z_0_yyzz_xyz[k] * cd_y[k] + g_z_0_yyzz_xyyz[k];

                g_z_0_yyyzz_xzz[k] = -g_z_0_yyzz_xzz[k] * cd_y[k] + g_z_0_yyzz_xyzz[k];

                g_z_0_yyyzz_yyy[k] = -g_z_0_yyzz_yyy[k] * cd_y[k] + g_z_0_yyzz_yyyy[k];

                g_z_0_yyyzz_yyz[k] = -g_z_0_yyzz_yyz[k] * cd_y[k] + g_z_0_yyzz_yyyz[k];

                g_z_0_yyyzz_yzz[k] = -g_z_0_yyzz_yzz[k] * cd_y[k] + g_z_0_yyzz_yyzz[k];

                g_z_0_yyyzz_zzz[k] = -g_z_0_yyzz_zzz[k] * cd_y[k] + g_z_0_yyzz_yzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 180);

            auto g_z_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 181);

            auto g_z_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 182);

            auto g_z_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 183);

            auto g_z_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 184);

            auto g_z_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 185);

            auto g_z_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 186);

            auto g_z_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 187);

            auto g_z_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 188);

            auto g_z_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 189);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzz_xxx, g_z_0_yyzzz_xxy, g_z_0_yyzzz_xxz, g_z_0_yyzzz_xyy, g_z_0_yyzzz_xyz, g_z_0_yyzzz_xzz, g_z_0_yyzzz_yyy, g_z_0_yyzzz_yyz, g_z_0_yyzzz_yzz, g_z_0_yyzzz_zzz, g_z_0_yzzz_xxx, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxy, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxz, g_z_0_yzzz_xyy, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xzz, g_z_0_yzzz_yyy, g_z_0_yzzz_yyyy, g_z_0_yzzz_yyyz, g_z_0_yzzz_yyz, g_z_0_yzzz_yyzz, g_z_0_yzzz_yzz, g_z_0_yzzz_yzzz, g_z_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_xxx[k] = -g_z_0_yzzz_xxx[k] * cd_y[k] + g_z_0_yzzz_xxxy[k];

                g_z_0_yyzzz_xxy[k] = -g_z_0_yzzz_xxy[k] * cd_y[k] + g_z_0_yzzz_xxyy[k];

                g_z_0_yyzzz_xxz[k] = -g_z_0_yzzz_xxz[k] * cd_y[k] + g_z_0_yzzz_xxyz[k];

                g_z_0_yyzzz_xyy[k] = -g_z_0_yzzz_xyy[k] * cd_y[k] + g_z_0_yzzz_xyyy[k];

                g_z_0_yyzzz_xyz[k] = -g_z_0_yzzz_xyz[k] * cd_y[k] + g_z_0_yzzz_xyyz[k];

                g_z_0_yyzzz_xzz[k] = -g_z_0_yzzz_xzz[k] * cd_y[k] + g_z_0_yzzz_xyzz[k];

                g_z_0_yyzzz_yyy[k] = -g_z_0_yzzz_yyy[k] * cd_y[k] + g_z_0_yzzz_yyyy[k];

                g_z_0_yyzzz_yyz[k] = -g_z_0_yzzz_yyz[k] * cd_y[k] + g_z_0_yzzz_yyyz[k];

                g_z_0_yyzzz_yzz[k] = -g_z_0_yzzz_yzz[k] * cd_y[k] + g_z_0_yzzz_yyzz[k];

                g_z_0_yyzzz_zzz[k] = -g_z_0_yzzz_zzz[k] * cd_y[k] + g_z_0_yzzz_yzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 190);

            auto g_z_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 191);

            auto g_z_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 192);

            auto g_z_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 193);

            auto g_z_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 194);

            auto g_z_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 195);

            auto g_z_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 196);

            auto g_z_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 197);

            auto g_z_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 198);

            auto g_z_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 199);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzz_xxx, g_z_0_yzzzz_xxy, g_z_0_yzzzz_xxz, g_z_0_yzzzz_xyy, g_z_0_yzzzz_xyz, g_z_0_yzzzz_xzz, g_z_0_yzzzz_yyy, g_z_0_yzzzz_yyz, g_z_0_yzzzz_yzz, g_z_0_yzzzz_zzz, g_z_0_zzzz_xxx, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxy, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxz, g_z_0_zzzz_xyy, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzz, g_z_0_zzzz_yyy, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_xxx[k] = -g_z_0_zzzz_xxx[k] * cd_y[k] + g_z_0_zzzz_xxxy[k];

                g_z_0_yzzzz_xxy[k] = -g_z_0_zzzz_xxy[k] * cd_y[k] + g_z_0_zzzz_xxyy[k];

                g_z_0_yzzzz_xxz[k] = -g_z_0_zzzz_xxz[k] * cd_y[k] + g_z_0_zzzz_xxyz[k];

                g_z_0_yzzzz_xyy[k] = -g_z_0_zzzz_xyy[k] * cd_y[k] + g_z_0_zzzz_xyyy[k];

                g_z_0_yzzzz_xyz[k] = -g_z_0_zzzz_xyz[k] * cd_y[k] + g_z_0_zzzz_xyyz[k];

                g_z_0_yzzzz_xzz[k] = -g_z_0_zzzz_xzz[k] * cd_y[k] + g_z_0_zzzz_xyzz[k];

                g_z_0_yzzzz_yyy[k] = -g_z_0_zzzz_yyy[k] * cd_y[k] + g_z_0_zzzz_yyyy[k];

                g_z_0_yzzzz_yyz[k] = -g_z_0_zzzz_yyz[k] * cd_y[k] + g_z_0_zzzz_yyyz[k];

                g_z_0_yzzzz_yzz[k] = -g_z_0_zzzz_yzz[k] * cd_y[k] + g_z_0_zzzz_yyzz[k];

                g_z_0_yzzzz_zzz[k] = -g_z_0_zzzz_zzz[k] * cd_y[k] + g_z_0_zzzz_yzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 200);

            auto g_z_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 201);

            auto g_z_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 202);

            auto g_z_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 203);

            auto g_z_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 204);

            auto g_z_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 205);

            auto g_z_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 206);

            auto g_z_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 207);

            auto g_z_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 208);

            auto g_z_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_z, g_z_0_zzzz_xxx, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzz, g_z_0_zzzz_zzzz, g_z_0_zzzzz_xxx, g_z_0_zzzzz_xxy, g_z_0_zzzzz_xxz, g_z_0_zzzzz_xyy, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xzz, g_z_0_zzzzz_yyy, g_z_0_zzzzz_yyz, g_z_0_zzzzz_yzz, g_z_0_zzzzz_zzz, g_zzzz_xxx, g_zzzz_xxy, g_zzzz_xxz, g_zzzz_xyy, g_zzzz_xyz, g_zzzz_xzz, g_zzzz_yyy, g_zzzz_yyz, g_zzzz_yzz, g_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_xxx[k] = -g_zzzz_xxx[k] - g_z_0_zzzz_xxx[k] * cd_z[k] + g_z_0_zzzz_xxxz[k];

                g_z_0_zzzzz_xxy[k] = -g_zzzz_xxy[k] - g_z_0_zzzz_xxy[k] * cd_z[k] + g_z_0_zzzz_xxyz[k];

                g_z_0_zzzzz_xxz[k] = -g_zzzz_xxz[k] - g_z_0_zzzz_xxz[k] * cd_z[k] + g_z_0_zzzz_xxzz[k];

                g_z_0_zzzzz_xyy[k] = -g_zzzz_xyy[k] - g_z_0_zzzz_xyy[k] * cd_z[k] + g_z_0_zzzz_xyyz[k];

                g_z_0_zzzzz_xyz[k] = -g_zzzz_xyz[k] - g_z_0_zzzz_xyz[k] * cd_z[k] + g_z_0_zzzz_xyzz[k];

                g_z_0_zzzzz_xzz[k] = -g_zzzz_xzz[k] - g_z_0_zzzz_xzz[k] * cd_z[k] + g_z_0_zzzz_xzzz[k];

                g_z_0_zzzzz_yyy[k] = -g_zzzz_yyy[k] - g_z_0_zzzz_yyy[k] * cd_z[k] + g_z_0_zzzz_yyyz[k];

                g_z_0_zzzzz_yyz[k] = -g_zzzz_yyz[k] - g_z_0_zzzz_yyz[k] * cd_z[k] + g_z_0_zzzz_yyzz[k];

                g_z_0_zzzzz_yzz[k] = -g_zzzz_yzz[k] - g_z_0_zzzz_yzz[k] * cd_z[k] + g_z_0_zzzz_yzzz[k];

                g_z_0_zzzzz_zzz[k] = -g_zzzz_zzz[k] - g_z_0_zzzz_zzz[k] * cd_z[k] + g_z_0_zzzz_zzzz[k];
            }
        }
    }
}

} // erirec namespace

