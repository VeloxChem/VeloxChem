#include "ElectronRepulsionGeom0010ContrRecXXHG.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxhg(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxhg,
                                            const size_t idx_xxgg,
                                            const size_t idx_geom_10_xxgg,
                                            const size_t idx_geom_10_xxgh,
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
            /// Set up components of auxilary buffer : SSGG

            const auto gg_off = idx_xxgg + (i * bcomps + j) * 225;

            auto g_xxxx_xxxx = cbuffer.data(gg_off + 0);

            auto g_xxxx_xxxy = cbuffer.data(gg_off + 1);

            auto g_xxxx_xxxz = cbuffer.data(gg_off + 2);

            auto g_xxxx_xxyy = cbuffer.data(gg_off + 3);

            auto g_xxxx_xxyz = cbuffer.data(gg_off + 4);

            auto g_xxxx_xxzz = cbuffer.data(gg_off + 5);

            auto g_xxxx_xyyy = cbuffer.data(gg_off + 6);

            auto g_xxxx_xyyz = cbuffer.data(gg_off + 7);

            auto g_xxxx_xyzz = cbuffer.data(gg_off + 8);

            auto g_xxxx_xzzz = cbuffer.data(gg_off + 9);

            auto g_xxxx_yyyy = cbuffer.data(gg_off + 10);

            auto g_xxxx_yyyz = cbuffer.data(gg_off + 11);

            auto g_xxxx_yyzz = cbuffer.data(gg_off + 12);

            auto g_xxxx_yzzz = cbuffer.data(gg_off + 13);

            auto g_xxxx_zzzz = cbuffer.data(gg_off + 14);

            auto g_yyyy_xxxx = cbuffer.data(gg_off + 150);

            auto g_yyyy_xxxy = cbuffer.data(gg_off + 151);

            auto g_yyyy_xxxz = cbuffer.data(gg_off + 152);

            auto g_yyyy_xxyy = cbuffer.data(gg_off + 153);

            auto g_yyyy_xxyz = cbuffer.data(gg_off + 154);

            auto g_yyyy_xxzz = cbuffer.data(gg_off + 155);

            auto g_yyyy_xyyy = cbuffer.data(gg_off + 156);

            auto g_yyyy_xyyz = cbuffer.data(gg_off + 157);

            auto g_yyyy_xyzz = cbuffer.data(gg_off + 158);

            auto g_yyyy_xzzz = cbuffer.data(gg_off + 159);

            auto g_yyyy_yyyy = cbuffer.data(gg_off + 160);

            auto g_yyyy_yyyz = cbuffer.data(gg_off + 161);

            auto g_yyyy_yyzz = cbuffer.data(gg_off + 162);

            auto g_yyyy_yzzz = cbuffer.data(gg_off + 163);

            auto g_yyyy_zzzz = cbuffer.data(gg_off + 164);

            auto g_zzzz_xxxx = cbuffer.data(gg_off + 210);

            auto g_zzzz_xxxy = cbuffer.data(gg_off + 211);

            auto g_zzzz_xxxz = cbuffer.data(gg_off + 212);

            auto g_zzzz_xxyy = cbuffer.data(gg_off + 213);

            auto g_zzzz_xxyz = cbuffer.data(gg_off + 214);

            auto g_zzzz_xxzz = cbuffer.data(gg_off + 215);

            auto g_zzzz_xyyy = cbuffer.data(gg_off + 216);

            auto g_zzzz_xyyz = cbuffer.data(gg_off + 217);

            auto g_zzzz_xyzz = cbuffer.data(gg_off + 218);

            auto g_zzzz_xzzz = cbuffer.data(gg_off + 219);

            auto g_zzzz_yyyy = cbuffer.data(gg_off + 220);

            auto g_zzzz_yyyz = cbuffer.data(gg_off + 221);

            auto g_zzzz_yyzz = cbuffer.data(gg_off + 222);

            auto g_zzzz_yzzz = cbuffer.data(gg_off + 223);

            auto g_zzzz_zzzz = cbuffer.data(gg_off + 224);

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

            auto g_x_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 30);

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

            auto g_x_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 75);

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

            auto g_x_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 135);

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

            auto g_x_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps * bcomps + 210);

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

            auto g_y_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 10);

            auto g_y_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 11);

            auto g_y_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 12);

            auto g_y_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 13);

            auto g_y_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 14);

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

            auto g_y_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 25);

            auto g_y_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 26);

            auto g_y_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 27);

            auto g_y_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 28);

            auto g_y_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 29);

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

            auto g_y_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 40);

            auto g_y_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 41);

            auto g_y_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 42);

            auto g_y_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 43);

            auto g_y_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 44);

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

            auto g_y_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 55);

            auto g_y_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 56);

            auto g_y_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 57);

            auto g_y_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 58);

            auto g_y_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 59);

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

            auto g_y_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 70);

            auto g_y_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 71);

            auto g_y_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 72);

            auto g_y_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 73);

            auto g_y_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 74);

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

            auto g_y_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 85);

            auto g_y_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 86);

            auto g_y_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 87);

            auto g_y_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 88);

            auto g_y_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 89);

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

            auto g_y_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 100);

            auto g_y_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 101);

            auto g_y_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 102);

            auto g_y_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 103);

            auto g_y_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 104);

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

            auto g_y_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 115);

            auto g_y_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 116);

            auto g_y_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 117);

            auto g_y_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 118);

            auto g_y_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 119);

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

            auto g_y_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 130);

            auto g_y_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 131);

            auto g_y_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 132);

            auto g_y_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 133);

            auto g_y_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 134);

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

            auto g_y_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 145);

            auto g_y_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 146);

            auto g_y_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 147);

            auto g_y_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 148);

            auto g_y_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 149);

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

            auto g_y_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 175);

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

            auto g_y_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 190);

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

            auto g_y_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 205);

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

            auto g_y_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps * bcomps + 220);

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

            auto g_z_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 10);

            auto g_z_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 11);

            auto g_z_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 12);

            auto g_z_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 13);

            auto g_z_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 14);

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

            auto g_z_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 25);

            auto g_z_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 26);

            auto g_z_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 27);

            auto g_z_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 28);

            auto g_z_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 29);

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

            auto g_z_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 40);

            auto g_z_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 41);

            auto g_z_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 42);

            auto g_z_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 43);

            auto g_z_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 44);

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

            auto g_z_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 55);

            auto g_z_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 56);

            auto g_z_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 57);

            auto g_z_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 58);

            auto g_z_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 59);

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

            auto g_z_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 70);

            auto g_z_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 71);

            auto g_z_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 72);

            auto g_z_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 73);

            auto g_z_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 74);

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

            auto g_z_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 85);

            auto g_z_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 86);

            auto g_z_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 87);

            auto g_z_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 88);

            auto g_z_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 89);

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

            auto g_z_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 100);

            auto g_z_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 101);

            auto g_z_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 102);

            auto g_z_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 103);

            auto g_z_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 104);

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

            auto g_z_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 115);

            auto g_z_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 116);

            auto g_z_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 117);

            auto g_z_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 118);

            auto g_z_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 119);

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

            auto g_z_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 130);

            auto g_z_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 131);

            auto g_z_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 132);

            auto g_z_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 133);

            auto g_z_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 134);

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

            auto g_z_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 145);

            auto g_z_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 146);

            auto g_z_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 147);

            auto g_z_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 148);

            auto g_z_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 149);

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

            auto g_z_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 164);

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

            auto g_z_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 179);

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

            auto g_z_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 194);

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

            auto g_z_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps * bcomps + 209);

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

            /// Set up components of auxilary buffer : SSGH

            const auto gh_geom_10_off = idx_geom_10_xxgh + (i * bcomps + j) * 315;

            auto g_x_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxx_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxx_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxx_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxx_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxx_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxx_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xzzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xzzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xzzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xzzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xzzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xzzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_y_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 2);

            auto g_y_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 3);

            auto g_y_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 4);

            auto g_y_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 5);

            auto g_y_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 6);

            auto g_y_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 7);

            auto g_y_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 8);

            auto g_y_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 9);

            auto g_y_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 10);

            auto g_y_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 11);

            auto g_y_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 12);

            auto g_y_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 13);

            auto g_y_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 14);

            auto g_y_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 21);

            auto g_y_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 22);

            auto g_y_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 23);

            auto g_y_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 24);

            auto g_y_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 25);

            auto g_y_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 26);

            auto g_y_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 27);

            auto g_y_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 28);

            auto g_y_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 29);

            auto g_y_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 30);

            auto g_y_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 31);

            auto g_y_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 32);

            auto g_y_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 33);

            auto g_y_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 34);

            auto g_y_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 35);

            auto g_y_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 42);

            auto g_y_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 43);

            auto g_y_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 44);

            auto g_y_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 45);

            auto g_y_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 46);

            auto g_y_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 47);

            auto g_y_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 48);

            auto g_y_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 49);

            auto g_y_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 50);

            auto g_y_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 51);

            auto g_y_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 52);

            auto g_y_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 53);

            auto g_y_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 54);

            auto g_y_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 55);

            auto g_y_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 56);

            auto g_y_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 63);

            auto g_y_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 64);

            auto g_y_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 65);

            auto g_y_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 66);

            auto g_y_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 67);

            auto g_y_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 68);

            auto g_y_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 69);

            auto g_y_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 70);

            auto g_y_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 71);

            auto g_y_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 72);

            auto g_y_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 73);

            auto g_y_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 74);

            auto g_y_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 75);

            auto g_y_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 76);

            auto g_y_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 77);

            auto g_y_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 84);

            auto g_y_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 85);

            auto g_y_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 86);

            auto g_y_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 87);

            auto g_y_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 88);

            auto g_y_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 89);

            auto g_y_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 90);

            auto g_y_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 91);

            auto g_y_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 92);

            auto g_y_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 93);

            auto g_y_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 94);

            auto g_y_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 95);

            auto g_y_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 96);

            auto g_y_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 97);

            auto g_y_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 98);

            auto g_y_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 105);

            auto g_y_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 106);

            auto g_y_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 107);

            auto g_y_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 108);

            auto g_y_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 109);

            auto g_y_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 110);

            auto g_y_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 111);

            auto g_y_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 112);

            auto g_y_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 113);

            auto g_y_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 114);

            auto g_y_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 115);

            auto g_y_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 116);

            auto g_y_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 117);

            auto g_y_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 118);

            auto g_y_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 119);

            auto g_y_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 126);

            auto g_y_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 127);

            auto g_y_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 128);

            auto g_y_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 129);

            auto g_y_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 130);

            auto g_y_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 131);

            auto g_y_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 132);

            auto g_y_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 133);

            auto g_y_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 134);

            auto g_y_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 135);

            auto g_y_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 136);

            auto g_y_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 137);

            auto g_y_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 138);

            auto g_y_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 139);

            auto g_y_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 140);

            auto g_y_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 147);

            auto g_y_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 148);

            auto g_y_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 149);

            auto g_y_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 150);

            auto g_y_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 151);

            auto g_y_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 152);

            auto g_y_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 153);

            auto g_y_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 154);

            auto g_y_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 155);

            auto g_y_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 156);

            auto g_y_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 157);

            auto g_y_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 158);

            auto g_y_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 159);

            auto g_y_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 160);

            auto g_y_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 161);

            auto g_y_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 168);

            auto g_y_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 169);

            auto g_y_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 170);

            auto g_y_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 171);

            auto g_y_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 172);

            auto g_y_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 173);

            auto g_y_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 174);

            auto g_y_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 175);

            auto g_y_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 176);

            auto g_y_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 177);

            auto g_y_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 178);

            auto g_y_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 179);

            auto g_y_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 180);

            auto g_y_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 181);

            auto g_y_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 182);

            auto g_y_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 189);

            auto g_y_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 190);

            auto g_y_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 191);

            auto g_y_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 192);

            auto g_y_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 193);

            auto g_y_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 194);

            auto g_y_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 195);

            auto g_y_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 196);

            auto g_y_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 197);

            auto g_y_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 198);

            auto g_y_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 199);

            auto g_y_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 200);

            auto g_y_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 201);

            auto g_y_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 202);

            auto g_y_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 203);

            auto g_y_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 210);

            auto g_y_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 211);

            auto g_y_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 212);

            auto g_y_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 213);

            auto g_y_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 214);

            auto g_y_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 215);

            auto g_y_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 216);

            auto g_y_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 217);

            auto g_y_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 218);

            auto g_y_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 219);

            auto g_y_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 220);

            auto g_y_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 221);

            auto g_y_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 222);

            auto g_y_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 223);

            auto g_y_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 224);

            auto g_y_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 225);

            auto g_y_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 226);

            auto g_y_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 227);

            auto g_y_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 228);

            auto g_y_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 229);

            auto g_y_0_yyyy_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 230);

            auto g_y_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 231);

            auto g_y_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 232);

            auto g_y_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 233);

            auto g_y_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 234);

            auto g_y_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 235);

            auto g_y_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 236);

            auto g_y_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 237);

            auto g_y_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 238);

            auto g_y_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 239);

            auto g_y_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 240);

            auto g_y_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 241);

            auto g_y_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 242);

            auto g_y_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 243);

            auto g_y_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 244);

            auto g_y_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 245);

            auto g_y_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 247);

            auto g_y_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 248);

            auto g_y_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 249);

            auto g_y_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 250);

            auto g_y_0_yyyz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 251);

            auto g_y_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 252);

            auto g_y_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 253);

            auto g_y_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 254);

            auto g_y_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 255);

            auto g_y_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 256);

            auto g_y_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 257);

            auto g_y_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 258);

            auto g_y_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 259);

            auto g_y_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 260);

            auto g_y_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 261);

            auto g_y_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 262);

            auto g_y_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 263);

            auto g_y_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 264);

            auto g_y_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 265);

            auto g_y_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 266);

            auto g_y_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 268);

            auto g_y_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 269);

            auto g_y_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 270);

            auto g_y_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 271);

            auto g_y_0_yyzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 272);

            auto g_y_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 273);

            auto g_y_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 274);

            auto g_y_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 275);

            auto g_y_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 276);

            auto g_y_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 277);

            auto g_y_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 278);

            auto g_y_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 279);

            auto g_y_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 280);

            auto g_y_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 281);

            auto g_y_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 282);

            auto g_y_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 283);

            auto g_y_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 284);

            auto g_y_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 285);

            auto g_y_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 286);

            auto g_y_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 287);

            auto g_y_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 289);

            auto g_y_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 290);

            auto g_y_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 291);

            auto g_y_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 292);

            auto g_y_0_yzzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 293);

            auto g_y_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 294);

            auto g_y_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 295);

            auto g_y_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 296);

            auto g_y_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 297);

            auto g_y_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 298);

            auto g_y_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 299);

            auto g_y_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 300);

            auto g_y_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 301);

            auto g_y_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 302);

            auto g_y_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 303);

            auto g_y_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 304);

            auto g_y_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 305);

            auto g_y_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 306);

            auto g_y_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 307);

            auto g_y_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 308);

            auto g_y_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 310);

            auto g_y_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 311);

            auto g_y_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 312);

            auto g_y_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 313);

            auto g_y_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 314);

            auto g_z_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 2);

            auto g_z_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 3);

            auto g_z_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 4);

            auto g_z_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 5);

            auto g_z_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 6);

            auto g_z_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 7);

            auto g_z_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 8);

            auto g_z_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 9);

            auto g_z_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 10);

            auto g_z_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 11);

            auto g_z_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 12);

            auto g_z_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 13);

            auto g_z_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 14);

            auto g_z_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 21);

            auto g_z_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 22);

            auto g_z_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 23);

            auto g_z_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 24);

            auto g_z_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 25);

            auto g_z_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 26);

            auto g_z_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 27);

            auto g_z_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 28);

            auto g_z_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 29);

            auto g_z_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 30);

            auto g_z_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 31);

            auto g_z_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 32);

            auto g_z_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 33);

            auto g_z_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 34);

            auto g_z_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 35);

            auto g_z_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 42);

            auto g_z_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 43);

            auto g_z_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 44);

            auto g_z_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 45);

            auto g_z_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 46);

            auto g_z_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 47);

            auto g_z_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 48);

            auto g_z_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 49);

            auto g_z_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 50);

            auto g_z_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 51);

            auto g_z_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 52);

            auto g_z_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 53);

            auto g_z_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 54);

            auto g_z_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 55);

            auto g_z_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 56);

            auto g_z_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 63);

            auto g_z_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 64);

            auto g_z_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 65);

            auto g_z_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 66);

            auto g_z_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 67);

            auto g_z_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 68);

            auto g_z_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 69);

            auto g_z_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 70);

            auto g_z_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 71);

            auto g_z_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 72);

            auto g_z_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 73);

            auto g_z_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 74);

            auto g_z_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 75);

            auto g_z_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 76);

            auto g_z_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 77);

            auto g_z_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 84);

            auto g_z_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 85);

            auto g_z_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 86);

            auto g_z_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 87);

            auto g_z_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 88);

            auto g_z_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 89);

            auto g_z_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 90);

            auto g_z_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 91);

            auto g_z_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 92);

            auto g_z_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 93);

            auto g_z_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 94);

            auto g_z_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 95);

            auto g_z_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 96);

            auto g_z_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 97);

            auto g_z_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 98);

            auto g_z_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 105);

            auto g_z_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 106);

            auto g_z_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 107);

            auto g_z_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 108);

            auto g_z_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 109);

            auto g_z_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 110);

            auto g_z_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 111);

            auto g_z_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 112);

            auto g_z_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 113);

            auto g_z_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 114);

            auto g_z_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 115);

            auto g_z_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 116);

            auto g_z_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 117);

            auto g_z_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 118);

            auto g_z_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 119);

            auto g_z_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 126);

            auto g_z_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 127);

            auto g_z_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 128);

            auto g_z_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 129);

            auto g_z_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 130);

            auto g_z_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 131);

            auto g_z_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 132);

            auto g_z_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 133);

            auto g_z_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 134);

            auto g_z_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 135);

            auto g_z_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 136);

            auto g_z_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 137);

            auto g_z_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 138);

            auto g_z_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 139);

            auto g_z_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 140);

            auto g_z_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 147);

            auto g_z_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 148);

            auto g_z_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 149);

            auto g_z_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 150);

            auto g_z_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 151);

            auto g_z_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 152);

            auto g_z_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 153);

            auto g_z_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 154);

            auto g_z_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 155);

            auto g_z_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 156);

            auto g_z_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 157);

            auto g_z_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 158);

            auto g_z_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 159);

            auto g_z_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 160);

            auto g_z_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 161);

            auto g_z_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 168);

            auto g_z_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 169);

            auto g_z_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 170);

            auto g_z_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 171);

            auto g_z_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 172);

            auto g_z_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 173);

            auto g_z_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 174);

            auto g_z_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 175);

            auto g_z_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 176);

            auto g_z_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 177);

            auto g_z_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 178);

            auto g_z_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 179);

            auto g_z_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 180);

            auto g_z_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 181);

            auto g_z_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 182);

            auto g_z_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 189);

            auto g_z_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 190);

            auto g_z_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 191);

            auto g_z_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 192);

            auto g_z_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 193);

            auto g_z_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 194);

            auto g_z_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 195);

            auto g_z_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 196);

            auto g_z_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 197);

            auto g_z_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 198);

            auto g_z_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 199);

            auto g_z_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 200);

            auto g_z_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 201);

            auto g_z_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 202);

            auto g_z_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 203);

            auto g_z_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 210);

            auto g_z_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 211);

            auto g_z_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 212);

            auto g_z_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 213);

            auto g_z_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 214);

            auto g_z_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 215);

            auto g_z_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 216);

            auto g_z_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 217);

            auto g_z_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 218);

            auto g_z_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 219);

            auto g_z_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 220);

            auto g_z_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 221);

            auto g_z_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 222);

            auto g_z_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 223);

            auto g_z_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 224);

            auto g_z_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 225);

            auto g_z_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 226);

            auto g_z_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 227);

            auto g_z_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 228);

            auto g_z_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 229);

            auto g_z_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 231);

            auto g_z_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 232);

            auto g_z_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 233);

            auto g_z_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 234);

            auto g_z_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 235);

            auto g_z_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 236);

            auto g_z_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 237);

            auto g_z_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 238);

            auto g_z_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 239);

            auto g_z_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 240);

            auto g_z_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 241);

            auto g_z_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 242);

            auto g_z_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 243);

            auto g_z_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 244);

            auto g_z_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 245);

            auto g_z_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 246);

            auto g_z_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 247);

            auto g_z_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 248);

            auto g_z_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 249);

            auto g_z_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 250);

            auto g_z_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 252);

            auto g_z_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 253);

            auto g_z_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 254);

            auto g_z_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 255);

            auto g_z_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 256);

            auto g_z_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 257);

            auto g_z_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 258);

            auto g_z_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 259);

            auto g_z_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 260);

            auto g_z_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 261);

            auto g_z_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 262);

            auto g_z_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 263);

            auto g_z_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 264);

            auto g_z_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 265);

            auto g_z_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 266);

            auto g_z_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 267);

            auto g_z_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 268);

            auto g_z_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 269);

            auto g_z_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 270);

            auto g_z_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 271);

            auto g_z_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 273);

            auto g_z_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 274);

            auto g_z_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 275);

            auto g_z_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 276);

            auto g_z_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 277);

            auto g_z_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 278);

            auto g_z_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 279);

            auto g_z_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 280);

            auto g_z_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 281);

            auto g_z_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 282);

            auto g_z_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 283);

            auto g_z_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 284);

            auto g_z_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 285);

            auto g_z_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 286);

            auto g_z_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 287);

            auto g_z_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 288);

            auto g_z_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 289);

            auto g_z_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 290);

            auto g_z_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 291);

            auto g_z_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 292);

            auto g_z_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 294);

            auto g_z_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 295);

            auto g_z_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 296);

            auto g_z_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 297);

            auto g_z_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 298);

            auto g_z_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 299);

            auto g_z_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 300);

            auto g_z_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 301);

            auto g_z_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 302);

            auto g_z_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 303);

            auto g_z_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 304);

            auto g_z_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 305);

            auto g_z_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 306);

            auto g_z_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 307);

            auto g_z_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 308);

            auto g_z_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 309);

            auto g_z_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 310);

            auto g_z_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 311);

            auto g_z_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 312);

            auto g_z_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 313);

            auto g_z_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 314);

            /// set up bra offset for contr_buffer_xxhg

            const auto hg_geom_10_off = idx_geom_10_xxhg + (i * bcomps + j) * 315;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxx_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxx_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxx_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxx_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxx_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxx_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxx_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxx_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxx_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxx_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxxx, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_zzzz, g_x_0_xxxxx_xxxx, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_yyyy, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_zzzz, g_xxxx_xxxx, g_xxxx_xxxy, g_xxxx_xxxz, g_xxxx_xxyy, g_xxxx_xxyz, g_xxxx_xxzz, g_xxxx_xyyy, g_xxxx_xyyz, g_xxxx_xyzz, g_xxxx_xzzz, g_xxxx_yyyy, g_xxxx_yyyz, g_xxxx_yyzz, g_xxxx_yzzz, g_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_xxxx[k] = -g_xxxx_xxxx[k] - g_x_0_xxxx_xxxx[k] * cd_x[k] + g_x_0_xxxx_xxxxx[k];

                g_x_0_xxxxx_xxxy[k] = -g_xxxx_xxxy[k] - g_x_0_xxxx_xxxy[k] * cd_x[k] + g_x_0_xxxx_xxxxy[k];

                g_x_0_xxxxx_xxxz[k] = -g_xxxx_xxxz[k] - g_x_0_xxxx_xxxz[k] * cd_x[k] + g_x_0_xxxx_xxxxz[k];

                g_x_0_xxxxx_xxyy[k] = -g_xxxx_xxyy[k] - g_x_0_xxxx_xxyy[k] * cd_x[k] + g_x_0_xxxx_xxxyy[k];

                g_x_0_xxxxx_xxyz[k] = -g_xxxx_xxyz[k] - g_x_0_xxxx_xxyz[k] * cd_x[k] + g_x_0_xxxx_xxxyz[k];

                g_x_0_xxxxx_xxzz[k] = -g_xxxx_xxzz[k] - g_x_0_xxxx_xxzz[k] * cd_x[k] + g_x_0_xxxx_xxxzz[k];

                g_x_0_xxxxx_xyyy[k] = -g_xxxx_xyyy[k] - g_x_0_xxxx_xyyy[k] * cd_x[k] + g_x_0_xxxx_xxyyy[k];

                g_x_0_xxxxx_xyyz[k] = -g_xxxx_xyyz[k] - g_x_0_xxxx_xyyz[k] * cd_x[k] + g_x_0_xxxx_xxyyz[k];

                g_x_0_xxxxx_xyzz[k] = -g_xxxx_xyzz[k] - g_x_0_xxxx_xyzz[k] * cd_x[k] + g_x_0_xxxx_xxyzz[k];

                g_x_0_xxxxx_xzzz[k] = -g_xxxx_xzzz[k] - g_x_0_xxxx_xzzz[k] * cd_x[k] + g_x_0_xxxx_xxzzz[k];

                g_x_0_xxxxx_yyyy[k] = -g_xxxx_yyyy[k] - g_x_0_xxxx_yyyy[k] * cd_x[k] + g_x_0_xxxx_xyyyy[k];

                g_x_0_xxxxx_yyyz[k] = -g_xxxx_yyyz[k] - g_x_0_xxxx_yyyz[k] * cd_x[k] + g_x_0_xxxx_xyyyz[k];

                g_x_0_xxxxx_yyzz[k] = -g_xxxx_yyzz[k] - g_x_0_xxxx_yyzz[k] * cd_x[k] + g_x_0_xxxx_xyyzz[k];

                g_x_0_xxxxx_yzzz[k] = -g_xxxx_yzzz[k] - g_x_0_xxxx_yzzz[k] * cd_x[k] + g_x_0_xxxx_xyzzz[k];

                g_x_0_xxxxx_zzzz[k] = -g_xxxx_zzzz[k] - g_x_0_xxxx_zzzz[k] * cd_x[k] + g_x_0_xxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxy_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxy_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxy_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxxy_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxy_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxxy_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxxy_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxy_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxy_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxy_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxy_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxy_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxy_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxy_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxy_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_y, g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyyy, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_zzzz, g_x_0_xxxxy_xxxx, g_x_0_xxxxy_xxxy, g_x_0_xxxxy_xxxz, g_x_0_xxxxy_xxyy, g_x_0_xxxxy_xxyz, g_x_0_xxxxy_xxzz, g_x_0_xxxxy_xyyy, g_x_0_xxxxy_xyyz, g_x_0_xxxxy_xyzz, g_x_0_xxxxy_xzzz, g_x_0_xxxxy_yyyy, g_x_0_xxxxy_yyyz, g_x_0_xxxxy_yyzz, g_x_0_xxxxy_yzzz, g_x_0_xxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_xxxx[k] = -g_x_0_xxxx_xxxx[k] * cd_y[k] + g_x_0_xxxx_xxxxy[k];

                g_x_0_xxxxy_xxxy[k] = -g_x_0_xxxx_xxxy[k] * cd_y[k] + g_x_0_xxxx_xxxyy[k];

                g_x_0_xxxxy_xxxz[k] = -g_x_0_xxxx_xxxz[k] * cd_y[k] + g_x_0_xxxx_xxxyz[k];

                g_x_0_xxxxy_xxyy[k] = -g_x_0_xxxx_xxyy[k] * cd_y[k] + g_x_0_xxxx_xxyyy[k];

                g_x_0_xxxxy_xxyz[k] = -g_x_0_xxxx_xxyz[k] * cd_y[k] + g_x_0_xxxx_xxyyz[k];

                g_x_0_xxxxy_xxzz[k] = -g_x_0_xxxx_xxzz[k] * cd_y[k] + g_x_0_xxxx_xxyzz[k];

                g_x_0_xxxxy_xyyy[k] = -g_x_0_xxxx_xyyy[k] * cd_y[k] + g_x_0_xxxx_xyyyy[k];

                g_x_0_xxxxy_xyyz[k] = -g_x_0_xxxx_xyyz[k] * cd_y[k] + g_x_0_xxxx_xyyyz[k];

                g_x_0_xxxxy_xyzz[k] = -g_x_0_xxxx_xyzz[k] * cd_y[k] + g_x_0_xxxx_xyyzz[k];

                g_x_0_xxxxy_xzzz[k] = -g_x_0_xxxx_xzzz[k] * cd_y[k] + g_x_0_xxxx_xyzzz[k];

                g_x_0_xxxxy_yyyy[k] = -g_x_0_xxxx_yyyy[k] * cd_y[k] + g_x_0_xxxx_yyyyy[k];

                g_x_0_xxxxy_yyyz[k] = -g_x_0_xxxx_yyyz[k] * cd_y[k] + g_x_0_xxxx_yyyyz[k];

                g_x_0_xxxxy_yyzz[k] = -g_x_0_xxxx_yyzz[k] * cd_y[k] + g_x_0_xxxx_yyyzz[k];

                g_x_0_xxxxy_yzzz[k] = -g_x_0_xxxx_yzzz[k] * cd_y[k] + g_x_0_xxxx_yyzzz[k];

                g_x_0_xxxxy_zzzz[k] = -g_x_0_xxxx_zzzz[k] * cd_y[k] + g_x_0_xxxx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxxz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxxz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxxz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxxz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxxz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxxz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxxz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxxz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxxz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxxz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxxz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxxz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxxz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxxz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_z, g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_zzzz, g_x_0_xxxx_zzzzz, g_x_0_xxxxz_xxxx, g_x_0_xxxxz_xxxy, g_x_0_xxxxz_xxxz, g_x_0_xxxxz_xxyy, g_x_0_xxxxz_xxyz, g_x_0_xxxxz_xxzz, g_x_0_xxxxz_xyyy, g_x_0_xxxxz_xyyz, g_x_0_xxxxz_xyzz, g_x_0_xxxxz_xzzz, g_x_0_xxxxz_yyyy, g_x_0_xxxxz_yyyz, g_x_0_xxxxz_yyzz, g_x_0_xxxxz_yzzz, g_x_0_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_xxxx[k] = -g_x_0_xxxx_xxxx[k] * cd_z[k] + g_x_0_xxxx_xxxxz[k];

                g_x_0_xxxxz_xxxy[k] = -g_x_0_xxxx_xxxy[k] * cd_z[k] + g_x_0_xxxx_xxxyz[k];

                g_x_0_xxxxz_xxxz[k] = -g_x_0_xxxx_xxxz[k] * cd_z[k] + g_x_0_xxxx_xxxzz[k];

                g_x_0_xxxxz_xxyy[k] = -g_x_0_xxxx_xxyy[k] * cd_z[k] + g_x_0_xxxx_xxyyz[k];

                g_x_0_xxxxz_xxyz[k] = -g_x_0_xxxx_xxyz[k] * cd_z[k] + g_x_0_xxxx_xxyzz[k];

                g_x_0_xxxxz_xxzz[k] = -g_x_0_xxxx_xxzz[k] * cd_z[k] + g_x_0_xxxx_xxzzz[k];

                g_x_0_xxxxz_xyyy[k] = -g_x_0_xxxx_xyyy[k] * cd_z[k] + g_x_0_xxxx_xyyyz[k];

                g_x_0_xxxxz_xyyz[k] = -g_x_0_xxxx_xyyz[k] * cd_z[k] + g_x_0_xxxx_xyyzz[k];

                g_x_0_xxxxz_xyzz[k] = -g_x_0_xxxx_xyzz[k] * cd_z[k] + g_x_0_xxxx_xyzzz[k];

                g_x_0_xxxxz_xzzz[k] = -g_x_0_xxxx_xzzz[k] * cd_z[k] + g_x_0_xxxx_xzzzz[k];

                g_x_0_xxxxz_yyyy[k] = -g_x_0_xxxx_yyyy[k] * cd_z[k] + g_x_0_xxxx_yyyyz[k];

                g_x_0_xxxxz_yyyz[k] = -g_x_0_xxxx_yyyz[k] * cd_z[k] + g_x_0_xxxx_yyyzz[k];

                g_x_0_xxxxz_yyzz[k] = -g_x_0_xxxx_yyzz[k] * cd_z[k] + g_x_0_xxxx_yyzzz[k];

                g_x_0_xxxxz_yzzz[k] = -g_x_0_xxxx_yzzz[k] * cd_z[k] + g_x_0_xxxx_yzzzz[k];

                g_x_0_xxxxz_zzzz[k] = -g_x_0_xxxx_zzzz[k] * cd_z[k] + g_x_0_xxxx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyy_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxyy_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxyy_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxyy_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxyy_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxyy_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxyy_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxyy_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxyy_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxyy_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxyy_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxyy_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxyy_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxyy_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxyy_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_y, g_x_0_xxxy_xxxx, g_x_0_xxxy_xxxxy, g_x_0_xxxy_xxxy, g_x_0_xxxy_xxxyy, g_x_0_xxxy_xxxyz, g_x_0_xxxy_xxxz, g_x_0_xxxy_xxyy, g_x_0_xxxy_xxyyy, g_x_0_xxxy_xxyyz, g_x_0_xxxy_xxyz, g_x_0_xxxy_xxyzz, g_x_0_xxxy_xxzz, g_x_0_xxxy_xyyy, g_x_0_xxxy_xyyyy, g_x_0_xxxy_xyyyz, g_x_0_xxxy_xyyz, g_x_0_xxxy_xyyzz, g_x_0_xxxy_xyzz, g_x_0_xxxy_xyzzz, g_x_0_xxxy_xzzz, g_x_0_xxxy_yyyy, g_x_0_xxxy_yyyyy, g_x_0_xxxy_yyyyz, g_x_0_xxxy_yyyz, g_x_0_xxxy_yyyzz, g_x_0_xxxy_yyzz, g_x_0_xxxy_yyzzz, g_x_0_xxxy_yzzz, g_x_0_xxxy_yzzzz, g_x_0_xxxy_zzzz, g_x_0_xxxyy_xxxx, g_x_0_xxxyy_xxxy, g_x_0_xxxyy_xxxz, g_x_0_xxxyy_xxyy, g_x_0_xxxyy_xxyz, g_x_0_xxxyy_xxzz, g_x_0_xxxyy_xyyy, g_x_0_xxxyy_xyyz, g_x_0_xxxyy_xyzz, g_x_0_xxxyy_xzzz, g_x_0_xxxyy_yyyy, g_x_0_xxxyy_yyyz, g_x_0_xxxyy_yyzz, g_x_0_xxxyy_yzzz, g_x_0_xxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_xxxx[k] = -g_x_0_xxxy_xxxx[k] * cd_y[k] + g_x_0_xxxy_xxxxy[k];

                g_x_0_xxxyy_xxxy[k] = -g_x_0_xxxy_xxxy[k] * cd_y[k] + g_x_0_xxxy_xxxyy[k];

                g_x_0_xxxyy_xxxz[k] = -g_x_0_xxxy_xxxz[k] * cd_y[k] + g_x_0_xxxy_xxxyz[k];

                g_x_0_xxxyy_xxyy[k] = -g_x_0_xxxy_xxyy[k] * cd_y[k] + g_x_0_xxxy_xxyyy[k];

                g_x_0_xxxyy_xxyz[k] = -g_x_0_xxxy_xxyz[k] * cd_y[k] + g_x_0_xxxy_xxyyz[k];

                g_x_0_xxxyy_xxzz[k] = -g_x_0_xxxy_xxzz[k] * cd_y[k] + g_x_0_xxxy_xxyzz[k];

                g_x_0_xxxyy_xyyy[k] = -g_x_0_xxxy_xyyy[k] * cd_y[k] + g_x_0_xxxy_xyyyy[k];

                g_x_0_xxxyy_xyyz[k] = -g_x_0_xxxy_xyyz[k] * cd_y[k] + g_x_0_xxxy_xyyyz[k];

                g_x_0_xxxyy_xyzz[k] = -g_x_0_xxxy_xyzz[k] * cd_y[k] + g_x_0_xxxy_xyyzz[k];

                g_x_0_xxxyy_xzzz[k] = -g_x_0_xxxy_xzzz[k] * cd_y[k] + g_x_0_xxxy_xyzzz[k];

                g_x_0_xxxyy_yyyy[k] = -g_x_0_xxxy_yyyy[k] * cd_y[k] + g_x_0_xxxy_yyyyy[k];

                g_x_0_xxxyy_yyyz[k] = -g_x_0_xxxy_yyyz[k] * cd_y[k] + g_x_0_xxxy_yyyyz[k];

                g_x_0_xxxyy_yyzz[k] = -g_x_0_xxxy_yyzz[k] * cd_y[k] + g_x_0_xxxy_yyyzz[k];

                g_x_0_xxxyy_yzzz[k] = -g_x_0_xxxy_yzzz[k] * cd_y[k] + g_x_0_xxxy_yyzzz[k];

                g_x_0_xxxyy_zzzz[k] = -g_x_0_xxxy_zzzz[k] * cd_y[k] + g_x_0_xxxy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxyz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxyz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxxyz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxyz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxyz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxyz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxyz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxyz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxyz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxyz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxyz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxyz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxyz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxyz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 74);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyz_xxxx, g_x_0_xxxyz_xxxy, g_x_0_xxxyz_xxxz, g_x_0_xxxyz_xxyy, g_x_0_xxxyz_xxyz, g_x_0_xxxyz_xxzz, g_x_0_xxxyz_xyyy, g_x_0_xxxyz_xyyz, g_x_0_xxxyz_xyzz, g_x_0_xxxyz_xzzz, g_x_0_xxxyz_yyyy, g_x_0_xxxyz_yyyz, g_x_0_xxxyz_yyzz, g_x_0_xxxyz_yzzz, g_x_0_xxxyz_zzzz, g_x_0_xxxz_xxxx, g_x_0_xxxz_xxxxy, g_x_0_xxxz_xxxy, g_x_0_xxxz_xxxyy, g_x_0_xxxz_xxxyz, g_x_0_xxxz_xxxz, g_x_0_xxxz_xxyy, g_x_0_xxxz_xxyyy, g_x_0_xxxz_xxyyz, g_x_0_xxxz_xxyz, g_x_0_xxxz_xxyzz, g_x_0_xxxz_xxzz, g_x_0_xxxz_xyyy, g_x_0_xxxz_xyyyy, g_x_0_xxxz_xyyyz, g_x_0_xxxz_xyyz, g_x_0_xxxz_xyyzz, g_x_0_xxxz_xyzz, g_x_0_xxxz_xyzzz, g_x_0_xxxz_xzzz, g_x_0_xxxz_yyyy, g_x_0_xxxz_yyyyy, g_x_0_xxxz_yyyyz, g_x_0_xxxz_yyyz, g_x_0_xxxz_yyyzz, g_x_0_xxxz_yyzz, g_x_0_xxxz_yyzzz, g_x_0_xxxz_yzzz, g_x_0_xxxz_yzzzz, g_x_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_xxxx[k] = -g_x_0_xxxz_xxxx[k] * cd_y[k] + g_x_0_xxxz_xxxxy[k];

                g_x_0_xxxyz_xxxy[k] = -g_x_0_xxxz_xxxy[k] * cd_y[k] + g_x_0_xxxz_xxxyy[k];

                g_x_0_xxxyz_xxxz[k] = -g_x_0_xxxz_xxxz[k] * cd_y[k] + g_x_0_xxxz_xxxyz[k];

                g_x_0_xxxyz_xxyy[k] = -g_x_0_xxxz_xxyy[k] * cd_y[k] + g_x_0_xxxz_xxyyy[k];

                g_x_0_xxxyz_xxyz[k] = -g_x_0_xxxz_xxyz[k] * cd_y[k] + g_x_0_xxxz_xxyyz[k];

                g_x_0_xxxyz_xxzz[k] = -g_x_0_xxxz_xxzz[k] * cd_y[k] + g_x_0_xxxz_xxyzz[k];

                g_x_0_xxxyz_xyyy[k] = -g_x_0_xxxz_xyyy[k] * cd_y[k] + g_x_0_xxxz_xyyyy[k];

                g_x_0_xxxyz_xyyz[k] = -g_x_0_xxxz_xyyz[k] * cd_y[k] + g_x_0_xxxz_xyyyz[k];

                g_x_0_xxxyz_xyzz[k] = -g_x_0_xxxz_xyzz[k] * cd_y[k] + g_x_0_xxxz_xyyzz[k];

                g_x_0_xxxyz_xzzz[k] = -g_x_0_xxxz_xzzz[k] * cd_y[k] + g_x_0_xxxz_xyzzz[k];

                g_x_0_xxxyz_yyyy[k] = -g_x_0_xxxz_yyyy[k] * cd_y[k] + g_x_0_xxxz_yyyyy[k];

                g_x_0_xxxyz_yyyz[k] = -g_x_0_xxxz_yyyz[k] * cd_y[k] + g_x_0_xxxz_yyyyz[k];

                g_x_0_xxxyz_yyzz[k] = -g_x_0_xxxz_yyzz[k] * cd_y[k] + g_x_0_xxxz_yyyzz[k];

                g_x_0_xxxyz_yzzz[k] = -g_x_0_xxxz_yzzz[k] * cd_y[k] + g_x_0_xxxz_yyzzz[k];

                g_x_0_xxxyz_zzzz[k] = -g_x_0_xxxz_zzzz[k] * cd_y[k] + g_x_0_xxxz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxxzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxxzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxxzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxxzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxxzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxxzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_z, g_x_0_xxxz_xxxx, g_x_0_xxxz_xxxxz, g_x_0_xxxz_xxxy, g_x_0_xxxz_xxxyz, g_x_0_xxxz_xxxz, g_x_0_xxxz_xxxzz, g_x_0_xxxz_xxyy, g_x_0_xxxz_xxyyz, g_x_0_xxxz_xxyz, g_x_0_xxxz_xxyzz, g_x_0_xxxz_xxzz, g_x_0_xxxz_xxzzz, g_x_0_xxxz_xyyy, g_x_0_xxxz_xyyyz, g_x_0_xxxz_xyyz, g_x_0_xxxz_xyyzz, g_x_0_xxxz_xyzz, g_x_0_xxxz_xyzzz, g_x_0_xxxz_xzzz, g_x_0_xxxz_xzzzz, g_x_0_xxxz_yyyy, g_x_0_xxxz_yyyyz, g_x_0_xxxz_yyyz, g_x_0_xxxz_yyyzz, g_x_0_xxxz_yyzz, g_x_0_xxxz_yyzzz, g_x_0_xxxz_yzzz, g_x_0_xxxz_yzzzz, g_x_0_xxxz_zzzz, g_x_0_xxxz_zzzzz, g_x_0_xxxzz_xxxx, g_x_0_xxxzz_xxxy, g_x_0_xxxzz_xxxz, g_x_0_xxxzz_xxyy, g_x_0_xxxzz_xxyz, g_x_0_xxxzz_xxzz, g_x_0_xxxzz_xyyy, g_x_0_xxxzz_xyyz, g_x_0_xxxzz_xyzz, g_x_0_xxxzz_xzzz, g_x_0_xxxzz_yyyy, g_x_0_xxxzz_yyyz, g_x_0_xxxzz_yyzz, g_x_0_xxxzz_yzzz, g_x_0_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_xxxx[k] = -g_x_0_xxxz_xxxx[k] * cd_z[k] + g_x_0_xxxz_xxxxz[k];

                g_x_0_xxxzz_xxxy[k] = -g_x_0_xxxz_xxxy[k] * cd_z[k] + g_x_0_xxxz_xxxyz[k];

                g_x_0_xxxzz_xxxz[k] = -g_x_0_xxxz_xxxz[k] * cd_z[k] + g_x_0_xxxz_xxxzz[k];

                g_x_0_xxxzz_xxyy[k] = -g_x_0_xxxz_xxyy[k] * cd_z[k] + g_x_0_xxxz_xxyyz[k];

                g_x_0_xxxzz_xxyz[k] = -g_x_0_xxxz_xxyz[k] * cd_z[k] + g_x_0_xxxz_xxyzz[k];

                g_x_0_xxxzz_xxzz[k] = -g_x_0_xxxz_xxzz[k] * cd_z[k] + g_x_0_xxxz_xxzzz[k];

                g_x_0_xxxzz_xyyy[k] = -g_x_0_xxxz_xyyy[k] * cd_z[k] + g_x_0_xxxz_xyyyz[k];

                g_x_0_xxxzz_xyyz[k] = -g_x_0_xxxz_xyyz[k] * cd_z[k] + g_x_0_xxxz_xyyzz[k];

                g_x_0_xxxzz_xyzz[k] = -g_x_0_xxxz_xyzz[k] * cd_z[k] + g_x_0_xxxz_xyzzz[k];

                g_x_0_xxxzz_xzzz[k] = -g_x_0_xxxz_xzzz[k] * cd_z[k] + g_x_0_xxxz_xzzzz[k];

                g_x_0_xxxzz_yyyy[k] = -g_x_0_xxxz_yyyy[k] * cd_z[k] + g_x_0_xxxz_yyyyz[k];

                g_x_0_xxxzz_yyyz[k] = -g_x_0_xxxz_yyyz[k] * cd_z[k] + g_x_0_xxxz_yyyzz[k];

                g_x_0_xxxzz_yyzz[k] = -g_x_0_xxxz_yyzz[k] * cd_z[k] + g_x_0_xxxz_yyzzz[k];

                g_x_0_xxxzz_yzzz[k] = -g_x_0_xxxz_yzzz[k] * cd_z[k] + g_x_0_xxxz_yzzzz[k];

                g_x_0_xxxzz_zzzz[k] = -g_x_0_xxxz_zzzz[k] * cd_z[k] + g_x_0_xxxz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyy_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxyyy_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxyyy_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxyyy_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxyyy_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxyyy_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxyyy_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxyyy_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxyyy_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxyyy_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxyyy_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxyyy_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxyyy_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxyyy_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxyyy_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_y, g_x_0_xxyy_xxxx, g_x_0_xxyy_xxxxy, g_x_0_xxyy_xxxy, g_x_0_xxyy_xxxyy, g_x_0_xxyy_xxxyz, g_x_0_xxyy_xxxz, g_x_0_xxyy_xxyy, g_x_0_xxyy_xxyyy, g_x_0_xxyy_xxyyz, g_x_0_xxyy_xxyz, g_x_0_xxyy_xxyzz, g_x_0_xxyy_xxzz, g_x_0_xxyy_xyyy, g_x_0_xxyy_xyyyy, g_x_0_xxyy_xyyyz, g_x_0_xxyy_xyyz, g_x_0_xxyy_xyyzz, g_x_0_xxyy_xyzz, g_x_0_xxyy_xyzzz, g_x_0_xxyy_xzzz, g_x_0_xxyy_yyyy, g_x_0_xxyy_yyyyy, g_x_0_xxyy_yyyyz, g_x_0_xxyy_yyyz, g_x_0_xxyy_yyyzz, g_x_0_xxyy_yyzz, g_x_0_xxyy_yyzzz, g_x_0_xxyy_yzzz, g_x_0_xxyy_yzzzz, g_x_0_xxyy_zzzz, g_x_0_xxyyy_xxxx, g_x_0_xxyyy_xxxy, g_x_0_xxyyy_xxxz, g_x_0_xxyyy_xxyy, g_x_0_xxyyy_xxyz, g_x_0_xxyyy_xxzz, g_x_0_xxyyy_xyyy, g_x_0_xxyyy_xyyz, g_x_0_xxyyy_xyzz, g_x_0_xxyyy_xzzz, g_x_0_xxyyy_yyyy, g_x_0_xxyyy_yyyz, g_x_0_xxyyy_yyzz, g_x_0_xxyyy_yzzz, g_x_0_xxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_xxxx[k] = -g_x_0_xxyy_xxxx[k] * cd_y[k] + g_x_0_xxyy_xxxxy[k];

                g_x_0_xxyyy_xxxy[k] = -g_x_0_xxyy_xxxy[k] * cd_y[k] + g_x_0_xxyy_xxxyy[k];

                g_x_0_xxyyy_xxxz[k] = -g_x_0_xxyy_xxxz[k] * cd_y[k] + g_x_0_xxyy_xxxyz[k];

                g_x_0_xxyyy_xxyy[k] = -g_x_0_xxyy_xxyy[k] * cd_y[k] + g_x_0_xxyy_xxyyy[k];

                g_x_0_xxyyy_xxyz[k] = -g_x_0_xxyy_xxyz[k] * cd_y[k] + g_x_0_xxyy_xxyyz[k];

                g_x_0_xxyyy_xxzz[k] = -g_x_0_xxyy_xxzz[k] * cd_y[k] + g_x_0_xxyy_xxyzz[k];

                g_x_0_xxyyy_xyyy[k] = -g_x_0_xxyy_xyyy[k] * cd_y[k] + g_x_0_xxyy_xyyyy[k];

                g_x_0_xxyyy_xyyz[k] = -g_x_0_xxyy_xyyz[k] * cd_y[k] + g_x_0_xxyy_xyyyz[k];

                g_x_0_xxyyy_xyzz[k] = -g_x_0_xxyy_xyzz[k] * cd_y[k] + g_x_0_xxyy_xyyzz[k];

                g_x_0_xxyyy_xzzz[k] = -g_x_0_xxyy_xzzz[k] * cd_y[k] + g_x_0_xxyy_xyzzz[k];

                g_x_0_xxyyy_yyyy[k] = -g_x_0_xxyy_yyyy[k] * cd_y[k] + g_x_0_xxyy_yyyyy[k];

                g_x_0_xxyyy_yyyz[k] = -g_x_0_xxyy_yyyz[k] * cd_y[k] + g_x_0_xxyy_yyyyz[k];

                g_x_0_xxyyy_yyzz[k] = -g_x_0_xxyy_yyzz[k] * cd_y[k] + g_x_0_xxyy_yyyzz[k];

                g_x_0_xxyyy_yzzz[k] = -g_x_0_xxyy_yzzz[k] * cd_y[k] + g_x_0_xxyy_yyzzz[k];

                g_x_0_xxyyy_zzzz[k] = -g_x_0_xxyy_zzzz[k] * cd_y[k] + g_x_0_xxyy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxyyz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxyyz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxyyz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxyyz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxyyz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxyyz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxyyz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxyyz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxyyz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxyyz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxyyz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxyyz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxyyz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxyyz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyz_xxxx, g_x_0_xxyyz_xxxy, g_x_0_xxyyz_xxxz, g_x_0_xxyyz_xxyy, g_x_0_xxyyz_xxyz, g_x_0_xxyyz_xxzz, g_x_0_xxyyz_xyyy, g_x_0_xxyyz_xyyz, g_x_0_xxyyz_xyzz, g_x_0_xxyyz_xzzz, g_x_0_xxyyz_yyyy, g_x_0_xxyyz_yyyz, g_x_0_xxyyz_yyzz, g_x_0_xxyyz_yzzz, g_x_0_xxyyz_zzzz, g_x_0_xxyz_xxxx, g_x_0_xxyz_xxxxy, g_x_0_xxyz_xxxy, g_x_0_xxyz_xxxyy, g_x_0_xxyz_xxxyz, g_x_0_xxyz_xxxz, g_x_0_xxyz_xxyy, g_x_0_xxyz_xxyyy, g_x_0_xxyz_xxyyz, g_x_0_xxyz_xxyz, g_x_0_xxyz_xxyzz, g_x_0_xxyz_xxzz, g_x_0_xxyz_xyyy, g_x_0_xxyz_xyyyy, g_x_0_xxyz_xyyyz, g_x_0_xxyz_xyyz, g_x_0_xxyz_xyyzz, g_x_0_xxyz_xyzz, g_x_0_xxyz_xyzzz, g_x_0_xxyz_xzzz, g_x_0_xxyz_yyyy, g_x_0_xxyz_yyyyy, g_x_0_xxyz_yyyyz, g_x_0_xxyz_yyyz, g_x_0_xxyz_yyyzz, g_x_0_xxyz_yyzz, g_x_0_xxyz_yyzzz, g_x_0_xxyz_yzzz, g_x_0_xxyz_yzzzz, g_x_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_xxxx[k] = -g_x_0_xxyz_xxxx[k] * cd_y[k] + g_x_0_xxyz_xxxxy[k];

                g_x_0_xxyyz_xxxy[k] = -g_x_0_xxyz_xxxy[k] * cd_y[k] + g_x_0_xxyz_xxxyy[k];

                g_x_0_xxyyz_xxxz[k] = -g_x_0_xxyz_xxxz[k] * cd_y[k] + g_x_0_xxyz_xxxyz[k];

                g_x_0_xxyyz_xxyy[k] = -g_x_0_xxyz_xxyy[k] * cd_y[k] + g_x_0_xxyz_xxyyy[k];

                g_x_0_xxyyz_xxyz[k] = -g_x_0_xxyz_xxyz[k] * cd_y[k] + g_x_0_xxyz_xxyyz[k];

                g_x_0_xxyyz_xxzz[k] = -g_x_0_xxyz_xxzz[k] * cd_y[k] + g_x_0_xxyz_xxyzz[k];

                g_x_0_xxyyz_xyyy[k] = -g_x_0_xxyz_xyyy[k] * cd_y[k] + g_x_0_xxyz_xyyyy[k];

                g_x_0_xxyyz_xyyz[k] = -g_x_0_xxyz_xyyz[k] * cd_y[k] + g_x_0_xxyz_xyyyz[k];

                g_x_0_xxyyz_xyzz[k] = -g_x_0_xxyz_xyzz[k] * cd_y[k] + g_x_0_xxyz_xyyzz[k];

                g_x_0_xxyyz_xzzz[k] = -g_x_0_xxyz_xzzz[k] * cd_y[k] + g_x_0_xxyz_xyzzz[k];

                g_x_0_xxyyz_yyyy[k] = -g_x_0_xxyz_yyyy[k] * cd_y[k] + g_x_0_xxyz_yyyyy[k];

                g_x_0_xxyyz_yyyz[k] = -g_x_0_xxyz_yyyz[k] * cd_y[k] + g_x_0_xxyz_yyyyz[k];

                g_x_0_xxyyz_yyzz[k] = -g_x_0_xxyz_yyzz[k] * cd_y[k] + g_x_0_xxyz_yyyzz[k];

                g_x_0_xxyyz_yzzz[k] = -g_x_0_xxyz_yzzz[k] * cd_y[k] + g_x_0_xxyz_yyzzz[k];

                g_x_0_xxyyz_zzzz[k] = -g_x_0_xxyz_zzzz[k] * cd_y[k] + g_x_0_xxyz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxyzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxyzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxyzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxyzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxyzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxyzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxyzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxyzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxyzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxyzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxyzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxyzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxyzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxyzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 134);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzz_xxxx, g_x_0_xxyzz_xxxy, g_x_0_xxyzz_xxxz, g_x_0_xxyzz_xxyy, g_x_0_xxyzz_xxyz, g_x_0_xxyzz_xxzz, g_x_0_xxyzz_xyyy, g_x_0_xxyzz_xyyz, g_x_0_xxyzz_xyzz, g_x_0_xxyzz_xzzz, g_x_0_xxyzz_yyyy, g_x_0_xxyzz_yyyz, g_x_0_xxyzz_yyzz, g_x_0_xxyzz_yzzz, g_x_0_xxyzz_zzzz, g_x_0_xxzz_xxxx, g_x_0_xxzz_xxxxy, g_x_0_xxzz_xxxy, g_x_0_xxzz_xxxyy, g_x_0_xxzz_xxxyz, g_x_0_xxzz_xxxz, g_x_0_xxzz_xxyy, g_x_0_xxzz_xxyyy, g_x_0_xxzz_xxyyz, g_x_0_xxzz_xxyz, g_x_0_xxzz_xxyzz, g_x_0_xxzz_xxzz, g_x_0_xxzz_xyyy, g_x_0_xxzz_xyyyy, g_x_0_xxzz_xyyyz, g_x_0_xxzz_xyyz, g_x_0_xxzz_xyyzz, g_x_0_xxzz_xyzz, g_x_0_xxzz_xyzzz, g_x_0_xxzz_xzzz, g_x_0_xxzz_yyyy, g_x_0_xxzz_yyyyy, g_x_0_xxzz_yyyyz, g_x_0_xxzz_yyyz, g_x_0_xxzz_yyyzz, g_x_0_xxzz_yyzz, g_x_0_xxzz_yyzzz, g_x_0_xxzz_yzzz, g_x_0_xxzz_yzzzz, g_x_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_xxxx[k] = -g_x_0_xxzz_xxxx[k] * cd_y[k] + g_x_0_xxzz_xxxxy[k];

                g_x_0_xxyzz_xxxy[k] = -g_x_0_xxzz_xxxy[k] * cd_y[k] + g_x_0_xxzz_xxxyy[k];

                g_x_0_xxyzz_xxxz[k] = -g_x_0_xxzz_xxxz[k] * cd_y[k] + g_x_0_xxzz_xxxyz[k];

                g_x_0_xxyzz_xxyy[k] = -g_x_0_xxzz_xxyy[k] * cd_y[k] + g_x_0_xxzz_xxyyy[k];

                g_x_0_xxyzz_xxyz[k] = -g_x_0_xxzz_xxyz[k] * cd_y[k] + g_x_0_xxzz_xxyyz[k];

                g_x_0_xxyzz_xxzz[k] = -g_x_0_xxzz_xxzz[k] * cd_y[k] + g_x_0_xxzz_xxyzz[k];

                g_x_0_xxyzz_xyyy[k] = -g_x_0_xxzz_xyyy[k] * cd_y[k] + g_x_0_xxzz_xyyyy[k];

                g_x_0_xxyzz_xyyz[k] = -g_x_0_xxzz_xyyz[k] * cd_y[k] + g_x_0_xxzz_xyyyz[k];

                g_x_0_xxyzz_xyzz[k] = -g_x_0_xxzz_xyzz[k] * cd_y[k] + g_x_0_xxzz_xyyzz[k];

                g_x_0_xxyzz_xzzz[k] = -g_x_0_xxzz_xzzz[k] * cd_y[k] + g_x_0_xxzz_xyzzz[k];

                g_x_0_xxyzz_yyyy[k] = -g_x_0_xxzz_yyyy[k] * cd_y[k] + g_x_0_xxzz_yyyyy[k];

                g_x_0_xxyzz_yyyz[k] = -g_x_0_xxzz_yyyz[k] * cd_y[k] + g_x_0_xxzz_yyyyz[k];

                g_x_0_xxyzz_yyzz[k] = -g_x_0_xxzz_yyzz[k] * cd_y[k] + g_x_0_xxzz_yyyzz[k];

                g_x_0_xxyzz_yzzz[k] = -g_x_0_xxzz_yzzz[k] * cd_y[k] + g_x_0_xxzz_yyzzz[k];

                g_x_0_xxyzz_zzzz[k] = -g_x_0_xxzz_zzzz[k] * cd_y[k] + g_x_0_xxzz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xxzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_z, g_x_0_xxzz_xxxx, g_x_0_xxzz_xxxxz, g_x_0_xxzz_xxxy, g_x_0_xxzz_xxxyz, g_x_0_xxzz_xxxz, g_x_0_xxzz_xxxzz, g_x_0_xxzz_xxyy, g_x_0_xxzz_xxyyz, g_x_0_xxzz_xxyz, g_x_0_xxzz_xxyzz, g_x_0_xxzz_xxzz, g_x_0_xxzz_xxzzz, g_x_0_xxzz_xyyy, g_x_0_xxzz_xyyyz, g_x_0_xxzz_xyyz, g_x_0_xxzz_xyyzz, g_x_0_xxzz_xyzz, g_x_0_xxzz_xyzzz, g_x_0_xxzz_xzzz, g_x_0_xxzz_xzzzz, g_x_0_xxzz_yyyy, g_x_0_xxzz_yyyyz, g_x_0_xxzz_yyyz, g_x_0_xxzz_yyyzz, g_x_0_xxzz_yyzz, g_x_0_xxzz_yyzzz, g_x_0_xxzz_yzzz, g_x_0_xxzz_yzzzz, g_x_0_xxzz_zzzz, g_x_0_xxzz_zzzzz, g_x_0_xxzzz_xxxx, g_x_0_xxzzz_xxxy, g_x_0_xxzzz_xxxz, g_x_0_xxzzz_xxyy, g_x_0_xxzzz_xxyz, g_x_0_xxzzz_xxzz, g_x_0_xxzzz_xyyy, g_x_0_xxzzz_xyyz, g_x_0_xxzzz_xyzz, g_x_0_xxzzz_xzzz, g_x_0_xxzzz_yyyy, g_x_0_xxzzz_yyyz, g_x_0_xxzzz_yyzz, g_x_0_xxzzz_yzzz, g_x_0_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_xxxx[k] = -g_x_0_xxzz_xxxx[k] * cd_z[k] + g_x_0_xxzz_xxxxz[k];

                g_x_0_xxzzz_xxxy[k] = -g_x_0_xxzz_xxxy[k] * cd_z[k] + g_x_0_xxzz_xxxyz[k];

                g_x_0_xxzzz_xxxz[k] = -g_x_0_xxzz_xxxz[k] * cd_z[k] + g_x_0_xxzz_xxxzz[k];

                g_x_0_xxzzz_xxyy[k] = -g_x_0_xxzz_xxyy[k] * cd_z[k] + g_x_0_xxzz_xxyyz[k];

                g_x_0_xxzzz_xxyz[k] = -g_x_0_xxzz_xxyz[k] * cd_z[k] + g_x_0_xxzz_xxyzz[k];

                g_x_0_xxzzz_xxzz[k] = -g_x_0_xxzz_xxzz[k] * cd_z[k] + g_x_0_xxzz_xxzzz[k];

                g_x_0_xxzzz_xyyy[k] = -g_x_0_xxzz_xyyy[k] * cd_z[k] + g_x_0_xxzz_xyyyz[k];

                g_x_0_xxzzz_xyyz[k] = -g_x_0_xxzz_xyyz[k] * cd_z[k] + g_x_0_xxzz_xyyzz[k];

                g_x_0_xxzzz_xyzz[k] = -g_x_0_xxzz_xyzz[k] * cd_z[k] + g_x_0_xxzz_xyzzz[k];

                g_x_0_xxzzz_xzzz[k] = -g_x_0_xxzz_xzzz[k] * cd_z[k] + g_x_0_xxzz_xzzzz[k];

                g_x_0_xxzzz_yyyy[k] = -g_x_0_xxzz_yyyy[k] * cd_z[k] + g_x_0_xxzz_yyyyz[k];

                g_x_0_xxzzz_yyyz[k] = -g_x_0_xxzz_yyyz[k] * cd_z[k] + g_x_0_xxzz_yyyzz[k];

                g_x_0_xxzzz_yyzz[k] = -g_x_0_xxzz_yyzz[k] * cd_z[k] + g_x_0_xxzz_yyzzz[k];

                g_x_0_xxzzz_yzzz[k] = -g_x_0_xxzz_yzzz[k] * cd_z[k] + g_x_0_xxzz_yzzzz[k];

                g_x_0_xxzzz_zzzz[k] = -g_x_0_xxzz_zzzz[k] * cd_z[k] + g_x_0_xxzz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyy_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xyyyy_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xyyyy_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xyyyy_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xyyyy_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xyyyy_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xyyyy_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xyyyy_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xyyyy_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xyyyy_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xyyyy_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xyyyy_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xyyyy_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xyyyy_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xyyyy_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 164);

            #pragma omp simd aligned(cd_y, g_x_0_xyyy_xxxx, g_x_0_xyyy_xxxxy, g_x_0_xyyy_xxxy, g_x_0_xyyy_xxxyy, g_x_0_xyyy_xxxyz, g_x_0_xyyy_xxxz, g_x_0_xyyy_xxyy, g_x_0_xyyy_xxyyy, g_x_0_xyyy_xxyyz, g_x_0_xyyy_xxyz, g_x_0_xyyy_xxyzz, g_x_0_xyyy_xxzz, g_x_0_xyyy_xyyy, g_x_0_xyyy_xyyyy, g_x_0_xyyy_xyyyz, g_x_0_xyyy_xyyz, g_x_0_xyyy_xyyzz, g_x_0_xyyy_xyzz, g_x_0_xyyy_xyzzz, g_x_0_xyyy_xzzz, g_x_0_xyyy_yyyy, g_x_0_xyyy_yyyyy, g_x_0_xyyy_yyyyz, g_x_0_xyyy_yyyz, g_x_0_xyyy_yyyzz, g_x_0_xyyy_yyzz, g_x_0_xyyy_yyzzz, g_x_0_xyyy_yzzz, g_x_0_xyyy_yzzzz, g_x_0_xyyy_zzzz, g_x_0_xyyyy_xxxx, g_x_0_xyyyy_xxxy, g_x_0_xyyyy_xxxz, g_x_0_xyyyy_xxyy, g_x_0_xyyyy_xxyz, g_x_0_xyyyy_xxzz, g_x_0_xyyyy_xyyy, g_x_0_xyyyy_xyyz, g_x_0_xyyyy_xyzz, g_x_0_xyyyy_xzzz, g_x_0_xyyyy_yyyy, g_x_0_xyyyy_yyyz, g_x_0_xyyyy_yyzz, g_x_0_xyyyy_yzzz, g_x_0_xyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_xxxx[k] = -g_x_0_xyyy_xxxx[k] * cd_y[k] + g_x_0_xyyy_xxxxy[k];

                g_x_0_xyyyy_xxxy[k] = -g_x_0_xyyy_xxxy[k] * cd_y[k] + g_x_0_xyyy_xxxyy[k];

                g_x_0_xyyyy_xxxz[k] = -g_x_0_xyyy_xxxz[k] * cd_y[k] + g_x_0_xyyy_xxxyz[k];

                g_x_0_xyyyy_xxyy[k] = -g_x_0_xyyy_xxyy[k] * cd_y[k] + g_x_0_xyyy_xxyyy[k];

                g_x_0_xyyyy_xxyz[k] = -g_x_0_xyyy_xxyz[k] * cd_y[k] + g_x_0_xyyy_xxyyz[k];

                g_x_0_xyyyy_xxzz[k] = -g_x_0_xyyy_xxzz[k] * cd_y[k] + g_x_0_xyyy_xxyzz[k];

                g_x_0_xyyyy_xyyy[k] = -g_x_0_xyyy_xyyy[k] * cd_y[k] + g_x_0_xyyy_xyyyy[k];

                g_x_0_xyyyy_xyyz[k] = -g_x_0_xyyy_xyyz[k] * cd_y[k] + g_x_0_xyyy_xyyyz[k];

                g_x_0_xyyyy_xyzz[k] = -g_x_0_xyyy_xyzz[k] * cd_y[k] + g_x_0_xyyy_xyyzz[k];

                g_x_0_xyyyy_xzzz[k] = -g_x_0_xyyy_xzzz[k] * cd_y[k] + g_x_0_xyyy_xyzzz[k];

                g_x_0_xyyyy_yyyy[k] = -g_x_0_xyyy_yyyy[k] * cd_y[k] + g_x_0_xyyy_yyyyy[k];

                g_x_0_xyyyy_yyyz[k] = -g_x_0_xyyy_yyyz[k] * cd_y[k] + g_x_0_xyyy_yyyyz[k];

                g_x_0_xyyyy_yyzz[k] = -g_x_0_xyyy_yyzz[k] * cd_y[k] + g_x_0_xyyy_yyyzz[k];

                g_x_0_xyyyy_yzzz[k] = -g_x_0_xyyy_yzzz[k] * cd_y[k] + g_x_0_xyyy_yyzzz[k];

                g_x_0_xyyyy_zzzz[k] = -g_x_0_xyyy_zzzz[k] * cd_y[k] + g_x_0_xyyy_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xyyyz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xyyyz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xyyyz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xyyyz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xyyyz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xyyyz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xyyyz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xyyyz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xyyyz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xyyyz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xyyyz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xyyyz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xyyyz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xyyyz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyz_xxxx, g_x_0_xyyyz_xxxy, g_x_0_xyyyz_xxxz, g_x_0_xyyyz_xxyy, g_x_0_xyyyz_xxyz, g_x_0_xyyyz_xxzz, g_x_0_xyyyz_xyyy, g_x_0_xyyyz_xyyz, g_x_0_xyyyz_xyzz, g_x_0_xyyyz_xzzz, g_x_0_xyyyz_yyyy, g_x_0_xyyyz_yyyz, g_x_0_xyyyz_yyzz, g_x_0_xyyyz_yzzz, g_x_0_xyyyz_zzzz, g_x_0_xyyz_xxxx, g_x_0_xyyz_xxxxy, g_x_0_xyyz_xxxy, g_x_0_xyyz_xxxyy, g_x_0_xyyz_xxxyz, g_x_0_xyyz_xxxz, g_x_0_xyyz_xxyy, g_x_0_xyyz_xxyyy, g_x_0_xyyz_xxyyz, g_x_0_xyyz_xxyz, g_x_0_xyyz_xxyzz, g_x_0_xyyz_xxzz, g_x_0_xyyz_xyyy, g_x_0_xyyz_xyyyy, g_x_0_xyyz_xyyyz, g_x_0_xyyz_xyyz, g_x_0_xyyz_xyyzz, g_x_0_xyyz_xyzz, g_x_0_xyyz_xyzzz, g_x_0_xyyz_xzzz, g_x_0_xyyz_yyyy, g_x_0_xyyz_yyyyy, g_x_0_xyyz_yyyyz, g_x_0_xyyz_yyyz, g_x_0_xyyz_yyyzz, g_x_0_xyyz_yyzz, g_x_0_xyyz_yyzzz, g_x_0_xyyz_yzzz, g_x_0_xyyz_yzzzz, g_x_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_xxxx[k] = -g_x_0_xyyz_xxxx[k] * cd_y[k] + g_x_0_xyyz_xxxxy[k];

                g_x_0_xyyyz_xxxy[k] = -g_x_0_xyyz_xxxy[k] * cd_y[k] + g_x_0_xyyz_xxxyy[k];

                g_x_0_xyyyz_xxxz[k] = -g_x_0_xyyz_xxxz[k] * cd_y[k] + g_x_0_xyyz_xxxyz[k];

                g_x_0_xyyyz_xxyy[k] = -g_x_0_xyyz_xxyy[k] * cd_y[k] + g_x_0_xyyz_xxyyy[k];

                g_x_0_xyyyz_xxyz[k] = -g_x_0_xyyz_xxyz[k] * cd_y[k] + g_x_0_xyyz_xxyyz[k];

                g_x_0_xyyyz_xxzz[k] = -g_x_0_xyyz_xxzz[k] * cd_y[k] + g_x_0_xyyz_xxyzz[k];

                g_x_0_xyyyz_xyyy[k] = -g_x_0_xyyz_xyyy[k] * cd_y[k] + g_x_0_xyyz_xyyyy[k];

                g_x_0_xyyyz_xyyz[k] = -g_x_0_xyyz_xyyz[k] * cd_y[k] + g_x_0_xyyz_xyyyz[k];

                g_x_0_xyyyz_xyzz[k] = -g_x_0_xyyz_xyzz[k] * cd_y[k] + g_x_0_xyyz_xyyzz[k];

                g_x_0_xyyyz_xzzz[k] = -g_x_0_xyyz_xzzz[k] * cd_y[k] + g_x_0_xyyz_xyzzz[k];

                g_x_0_xyyyz_yyyy[k] = -g_x_0_xyyz_yyyy[k] * cd_y[k] + g_x_0_xyyz_yyyyy[k];

                g_x_0_xyyyz_yyyz[k] = -g_x_0_xyyz_yyyz[k] * cd_y[k] + g_x_0_xyyz_yyyyz[k];

                g_x_0_xyyyz_yyzz[k] = -g_x_0_xyyz_yyzz[k] * cd_y[k] + g_x_0_xyyz_yyyzz[k];

                g_x_0_xyyyz_yzzz[k] = -g_x_0_xyyz_yzzz[k] * cd_y[k] + g_x_0_xyyz_yyzzz[k];

                g_x_0_xyyyz_zzzz[k] = -g_x_0_xyyz_zzzz[k] * cd_y[k] + g_x_0_xyyz_yzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xyyzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xyyzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xyyzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xyyzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xyyzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xyyzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xyyzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xyyzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xyyzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xyyzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xyyzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xyyzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xyyzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xyyzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 194);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzz_xxxx, g_x_0_xyyzz_xxxy, g_x_0_xyyzz_xxxz, g_x_0_xyyzz_xxyy, g_x_0_xyyzz_xxyz, g_x_0_xyyzz_xxzz, g_x_0_xyyzz_xyyy, g_x_0_xyyzz_xyyz, g_x_0_xyyzz_xyzz, g_x_0_xyyzz_xzzz, g_x_0_xyyzz_yyyy, g_x_0_xyyzz_yyyz, g_x_0_xyyzz_yyzz, g_x_0_xyyzz_yzzz, g_x_0_xyyzz_zzzz, g_x_0_xyzz_xxxx, g_x_0_xyzz_xxxxy, g_x_0_xyzz_xxxy, g_x_0_xyzz_xxxyy, g_x_0_xyzz_xxxyz, g_x_0_xyzz_xxxz, g_x_0_xyzz_xxyy, g_x_0_xyzz_xxyyy, g_x_0_xyzz_xxyyz, g_x_0_xyzz_xxyz, g_x_0_xyzz_xxyzz, g_x_0_xyzz_xxzz, g_x_0_xyzz_xyyy, g_x_0_xyzz_xyyyy, g_x_0_xyzz_xyyyz, g_x_0_xyzz_xyyz, g_x_0_xyzz_xyyzz, g_x_0_xyzz_xyzz, g_x_0_xyzz_xyzzz, g_x_0_xyzz_xzzz, g_x_0_xyzz_yyyy, g_x_0_xyzz_yyyyy, g_x_0_xyzz_yyyyz, g_x_0_xyzz_yyyz, g_x_0_xyzz_yyyzz, g_x_0_xyzz_yyzz, g_x_0_xyzz_yyzzz, g_x_0_xyzz_yzzz, g_x_0_xyzz_yzzzz, g_x_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_xxxx[k] = -g_x_0_xyzz_xxxx[k] * cd_y[k] + g_x_0_xyzz_xxxxy[k];

                g_x_0_xyyzz_xxxy[k] = -g_x_0_xyzz_xxxy[k] * cd_y[k] + g_x_0_xyzz_xxxyy[k];

                g_x_0_xyyzz_xxxz[k] = -g_x_0_xyzz_xxxz[k] * cd_y[k] + g_x_0_xyzz_xxxyz[k];

                g_x_0_xyyzz_xxyy[k] = -g_x_0_xyzz_xxyy[k] * cd_y[k] + g_x_0_xyzz_xxyyy[k];

                g_x_0_xyyzz_xxyz[k] = -g_x_0_xyzz_xxyz[k] * cd_y[k] + g_x_0_xyzz_xxyyz[k];

                g_x_0_xyyzz_xxzz[k] = -g_x_0_xyzz_xxzz[k] * cd_y[k] + g_x_0_xyzz_xxyzz[k];

                g_x_0_xyyzz_xyyy[k] = -g_x_0_xyzz_xyyy[k] * cd_y[k] + g_x_0_xyzz_xyyyy[k];

                g_x_0_xyyzz_xyyz[k] = -g_x_0_xyzz_xyyz[k] * cd_y[k] + g_x_0_xyzz_xyyyz[k];

                g_x_0_xyyzz_xyzz[k] = -g_x_0_xyzz_xyzz[k] * cd_y[k] + g_x_0_xyzz_xyyzz[k];

                g_x_0_xyyzz_xzzz[k] = -g_x_0_xyzz_xzzz[k] * cd_y[k] + g_x_0_xyzz_xyzzz[k];

                g_x_0_xyyzz_yyyy[k] = -g_x_0_xyzz_yyyy[k] * cd_y[k] + g_x_0_xyzz_yyyyy[k];

                g_x_0_xyyzz_yyyz[k] = -g_x_0_xyzz_yyyz[k] * cd_y[k] + g_x_0_xyzz_yyyyz[k];

                g_x_0_xyyzz_yyzz[k] = -g_x_0_xyzz_yyzz[k] * cd_y[k] + g_x_0_xyzz_yyyzz[k];

                g_x_0_xyyzz_yzzz[k] = -g_x_0_xyzz_yzzz[k] * cd_y[k] + g_x_0_xyzz_yyzzz[k];

                g_x_0_xyyzz_zzzz[k] = -g_x_0_xyzz_zzzz[k] * cd_y[k] + g_x_0_xyzz_yzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xyzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xyzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xyzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xyzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xyzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xyzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xyzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xyzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xyzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xyzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xyzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xyzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xyzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xyzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzz_xxxx, g_x_0_xyzzz_xxxy, g_x_0_xyzzz_xxxz, g_x_0_xyzzz_xxyy, g_x_0_xyzzz_xxyz, g_x_0_xyzzz_xxzz, g_x_0_xyzzz_xyyy, g_x_0_xyzzz_xyyz, g_x_0_xyzzz_xyzz, g_x_0_xyzzz_xzzz, g_x_0_xyzzz_yyyy, g_x_0_xyzzz_yyyz, g_x_0_xyzzz_yyzz, g_x_0_xyzzz_yzzz, g_x_0_xyzzz_zzzz, g_x_0_xzzz_xxxx, g_x_0_xzzz_xxxxy, g_x_0_xzzz_xxxy, g_x_0_xzzz_xxxyy, g_x_0_xzzz_xxxyz, g_x_0_xzzz_xxxz, g_x_0_xzzz_xxyy, g_x_0_xzzz_xxyyy, g_x_0_xzzz_xxyyz, g_x_0_xzzz_xxyz, g_x_0_xzzz_xxyzz, g_x_0_xzzz_xxzz, g_x_0_xzzz_xyyy, g_x_0_xzzz_xyyyy, g_x_0_xzzz_xyyyz, g_x_0_xzzz_xyyz, g_x_0_xzzz_xyyzz, g_x_0_xzzz_xyzz, g_x_0_xzzz_xyzzz, g_x_0_xzzz_xzzz, g_x_0_xzzz_yyyy, g_x_0_xzzz_yyyyy, g_x_0_xzzz_yyyyz, g_x_0_xzzz_yyyz, g_x_0_xzzz_yyyzz, g_x_0_xzzz_yyzz, g_x_0_xzzz_yyzzz, g_x_0_xzzz_yzzz, g_x_0_xzzz_yzzzz, g_x_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_xxxx[k] = -g_x_0_xzzz_xxxx[k] * cd_y[k] + g_x_0_xzzz_xxxxy[k];

                g_x_0_xyzzz_xxxy[k] = -g_x_0_xzzz_xxxy[k] * cd_y[k] + g_x_0_xzzz_xxxyy[k];

                g_x_0_xyzzz_xxxz[k] = -g_x_0_xzzz_xxxz[k] * cd_y[k] + g_x_0_xzzz_xxxyz[k];

                g_x_0_xyzzz_xxyy[k] = -g_x_0_xzzz_xxyy[k] * cd_y[k] + g_x_0_xzzz_xxyyy[k];

                g_x_0_xyzzz_xxyz[k] = -g_x_0_xzzz_xxyz[k] * cd_y[k] + g_x_0_xzzz_xxyyz[k];

                g_x_0_xyzzz_xxzz[k] = -g_x_0_xzzz_xxzz[k] * cd_y[k] + g_x_0_xzzz_xxyzz[k];

                g_x_0_xyzzz_xyyy[k] = -g_x_0_xzzz_xyyy[k] * cd_y[k] + g_x_0_xzzz_xyyyy[k];

                g_x_0_xyzzz_xyyz[k] = -g_x_0_xzzz_xyyz[k] * cd_y[k] + g_x_0_xzzz_xyyyz[k];

                g_x_0_xyzzz_xyzz[k] = -g_x_0_xzzz_xyzz[k] * cd_y[k] + g_x_0_xzzz_xyyzz[k];

                g_x_0_xyzzz_xzzz[k] = -g_x_0_xzzz_xzzz[k] * cd_y[k] + g_x_0_xzzz_xyzzz[k];

                g_x_0_xyzzz_yyyy[k] = -g_x_0_xzzz_yyyy[k] * cd_y[k] + g_x_0_xzzz_yyyyy[k];

                g_x_0_xyzzz_yyyz[k] = -g_x_0_xzzz_yyyz[k] * cd_y[k] + g_x_0_xzzz_yyyyz[k];

                g_x_0_xyzzz_yyzz[k] = -g_x_0_xzzz_yyzz[k] * cd_y[k] + g_x_0_xzzz_yyyzz[k];

                g_x_0_xyzzz_yzzz[k] = -g_x_0_xzzz_yzzz[k] * cd_y[k] + g_x_0_xzzz_yyzzz[k];

                g_x_0_xyzzz_zzzz[k] = -g_x_0_xzzz_zzzz[k] * cd_y[k] + g_x_0_xzzz_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xzzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xzzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xzzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xzzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xzzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xzzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xzzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xzzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xzzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xzzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xzzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xzzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xzzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_xzzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 224);

            #pragma omp simd aligned(cd_z, g_x_0_xzzz_xxxx, g_x_0_xzzz_xxxxz, g_x_0_xzzz_xxxy, g_x_0_xzzz_xxxyz, g_x_0_xzzz_xxxz, g_x_0_xzzz_xxxzz, g_x_0_xzzz_xxyy, g_x_0_xzzz_xxyyz, g_x_0_xzzz_xxyz, g_x_0_xzzz_xxyzz, g_x_0_xzzz_xxzz, g_x_0_xzzz_xxzzz, g_x_0_xzzz_xyyy, g_x_0_xzzz_xyyyz, g_x_0_xzzz_xyyz, g_x_0_xzzz_xyyzz, g_x_0_xzzz_xyzz, g_x_0_xzzz_xyzzz, g_x_0_xzzz_xzzz, g_x_0_xzzz_xzzzz, g_x_0_xzzz_yyyy, g_x_0_xzzz_yyyyz, g_x_0_xzzz_yyyz, g_x_0_xzzz_yyyzz, g_x_0_xzzz_yyzz, g_x_0_xzzz_yyzzz, g_x_0_xzzz_yzzz, g_x_0_xzzz_yzzzz, g_x_0_xzzz_zzzz, g_x_0_xzzz_zzzzz, g_x_0_xzzzz_xxxx, g_x_0_xzzzz_xxxy, g_x_0_xzzzz_xxxz, g_x_0_xzzzz_xxyy, g_x_0_xzzzz_xxyz, g_x_0_xzzzz_xxzz, g_x_0_xzzzz_xyyy, g_x_0_xzzzz_xyyz, g_x_0_xzzzz_xyzz, g_x_0_xzzzz_xzzz, g_x_0_xzzzz_yyyy, g_x_0_xzzzz_yyyz, g_x_0_xzzzz_yyzz, g_x_0_xzzzz_yzzz, g_x_0_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_xxxx[k] = -g_x_0_xzzz_xxxx[k] * cd_z[k] + g_x_0_xzzz_xxxxz[k];

                g_x_0_xzzzz_xxxy[k] = -g_x_0_xzzz_xxxy[k] * cd_z[k] + g_x_0_xzzz_xxxyz[k];

                g_x_0_xzzzz_xxxz[k] = -g_x_0_xzzz_xxxz[k] * cd_z[k] + g_x_0_xzzz_xxxzz[k];

                g_x_0_xzzzz_xxyy[k] = -g_x_0_xzzz_xxyy[k] * cd_z[k] + g_x_0_xzzz_xxyyz[k];

                g_x_0_xzzzz_xxyz[k] = -g_x_0_xzzz_xxyz[k] * cd_z[k] + g_x_0_xzzz_xxyzz[k];

                g_x_0_xzzzz_xxzz[k] = -g_x_0_xzzz_xxzz[k] * cd_z[k] + g_x_0_xzzz_xxzzz[k];

                g_x_0_xzzzz_xyyy[k] = -g_x_0_xzzz_xyyy[k] * cd_z[k] + g_x_0_xzzz_xyyyz[k];

                g_x_0_xzzzz_xyyz[k] = -g_x_0_xzzz_xyyz[k] * cd_z[k] + g_x_0_xzzz_xyyzz[k];

                g_x_0_xzzzz_xyzz[k] = -g_x_0_xzzz_xyzz[k] * cd_z[k] + g_x_0_xzzz_xyzzz[k];

                g_x_0_xzzzz_xzzz[k] = -g_x_0_xzzz_xzzz[k] * cd_z[k] + g_x_0_xzzz_xzzzz[k];

                g_x_0_xzzzz_yyyy[k] = -g_x_0_xzzz_yyyy[k] * cd_z[k] + g_x_0_xzzz_yyyyz[k];

                g_x_0_xzzzz_yyyz[k] = -g_x_0_xzzz_yyyz[k] * cd_z[k] + g_x_0_xzzz_yyyzz[k];

                g_x_0_xzzzz_yyzz[k] = -g_x_0_xzzz_yyzz[k] * cd_z[k] + g_x_0_xzzz_yyzzz[k];

                g_x_0_xzzzz_yzzz[k] = -g_x_0_xzzz_yzzz[k] * cd_z[k] + g_x_0_xzzz_yzzzz[k];

                g_x_0_xzzzz_zzzz[k] = -g_x_0_xzzz_zzzz[k] * cd_z[k] + g_x_0_xzzz_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyy_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_yyyyy_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_yyyyy_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_yyyyy_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_yyyyy_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_yyyyy_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_yyyyy_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_yyyyy_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_yyyyy_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_yyyyy_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_yyyyy_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_yyyyy_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_yyyyy_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_yyyyy_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_yyyyy_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 239);

            #pragma omp simd aligned(cd_y, g_x_0_yyyy_xxxx, g_x_0_yyyy_xxxxy, g_x_0_yyyy_xxxy, g_x_0_yyyy_xxxyy, g_x_0_yyyy_xxxyz, g_x_0_yyyy_xxxz, g_x_0_yyyy_xxyy, g_x_0_yyyy_xxyyy, g_x_0_yyyy_xxyyz, g_x_0_yyyy_xxyz, g_x_0_yyyy_xxyzz, g_x_0_yyyy_xxzz, g_x_0_yyyy_xyyy, g_x_0_yyyy_xyyyy, g_x_0_yyyy_xyyyz, g_x_0_yyyy_xyyz, g_x_0_yyyy_xyyzz, g_x_0_yyyy_xyzz, g_x_0_yyyy_xyzzz, g_x_0_yyyy_xzzz, g_x_0_yyyy_yyyy, g_x_0_yyyy_yyyyy, g_x_0_yyyy_yyyyz, g_x_0_yyyy_yyyz, g_x_0_yyyy_yyyzz, g_x_0_yyyy_yyzz, g_x_0_yyyy_yyzzz, g_x_0_yyyy_yzzz, g_x_0_yyyy_yzzzz, g_x_0_yyyy_zzzz, g_x_0_yyyyy_xxxx, g_x_0_yyyyy_xxxy, g_x_0_yyyyy_xxxz, g_x_0_yyyyy_xxyy, g_x_0_yyyyy_xxyz, g_x_0_yyyyy_xxzz, g_x_0_yyyyy_xyyy, g_x_0_yyyyy_xyyz, g_x_0_yyyyy_xyzz, g_x_0_yyyyy_xzzz, g_x_0_yyyyy_yyyy, g_x_0_yyyyy_yyyz, g_x_0_yyyyy_yyzz, g_x_0_yyyyy_yzzz, g_x_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_xxxx[k] = -g_x_0_yyyy_xxxx[k] * cd_y[k] + g_x_0_yyyy_xxxxy[k];

                g_x_0_yyyyy_xxxy[k] = -g_x_0_yyyy_xxxy[k] * cd_y[k] + g_x_0_yyyy_xxxyy[k];

                g_x_0_yyyyy_xxxz[k] = -g_x_0_yyyy_xxxz[k] * cd_y[k] + g_x_0_yyyy_xxxyz[k];

                g_x_0_yyyyy_xxyy[k] = -g_x_0_yyyy_xxyy[k] * cd_y[k] + g_x_0_yyyy_xxyyy[k];

                g_x_0_yyyyy_xxyz[k] = -g_x_0_yyyy_xxyz[k] * cd_y[k] + g_x_0_yyyy_xxyyz[k];

                g_x_0_yyyyy_xxzz[k] = -g_x_0_yyyy_xxzz[k] * cd_y[k] + g_x_0_yyyy_xxyzz[k];

                g_x_0_yyyyy_xyyy[k] = -g_x_0_yyyy_xyyy[k] * cd_y[k] + g_x_0_yyyy_xyyyy[k];

                g_x_0_yyyyy_xyyz[k] = -g_x_0_yyyy_xyyz[k] * cd_y[k] + g_x_0_yyyy_xyyyz[k];

                g_x_0_yyyyy_xyzz[k] = -g_x_0_yyyy_xyzz[k] * cd_y[k] + g_x_0_yyyy_xyyzz[k];

                g_x_0_yyyyy_xzzz[k] = -g_x_0_yyyy_xzzz[k] * cd_y[k] + g_x_0_yyyy_xyzzz[k];

                g_x_0_yyyyy_yyyy[k] = -g_x_0_yyyy_yyyy[k] * cd_y[k] + g_x_0_yyyy_yyyyy[k];

                g_x_0_yyyyy_yyyz[k] = -g_x_0_yyyy_yyyz[k] * cd_y[k] + g_x_0_yyyy_yyyyz[k];

                g_x_0_yyyyy_yyzz[k] = -g_x_0_yyyy_yyzz[k] * cd_y[k] + g_x_0_yyyy_yyyzz[k];

                g_x_0_yyyyy_yzzz[k] = -g_x_0_yyyy_yzzz[k] * cd_y[k] + g_x_0_yyyy_yyzzz[k];

                g_x_0_yyyyy_zzzz[k] = -g_x_0_yyyy_zzzz[k] * cd_y[k] + g_x_0_yyyy_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_yyyyz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_yyyyz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_yyyyz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_yyyyz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_yyyyz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_yyyyz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_yyyyz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_yyyyz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_yyyyz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_yyyyz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_yyyyz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_yyyyz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_yyyyz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_yyyyz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 254);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyz_xxxx, g_x_0_yyyyz_xxxy, g_x_0_yyyyz_xxxz, g_x_0_yyyyz_xxyy, g_x_0_yyyyz_xxyz, g_x_0_yyyyz_xxzz, g_x_0_yyyyz_xyyy, g_x_0_yyyyz_xyyz, g_x_0_yyyyz_xyzz, g_x_0_yyyyz_xzzz, g_x_0_yyyyz_yyyy, g_x_0_yyyyz_yyyz, g_x_0_yyyyz_yyzz, g_x_0_yyyyz_yzzz, g_x_0_yyyyz_zzzz, g_x_0_yyyz_xxxx, g_x_0_yyyz_xxxxy, g_x_0_yyyz_xxxy, g_x_0_yyyz_xxxyy, g_x_0_yyyz_xxxyz, g_x_0_yyyz_xxxz, g_x_0_yyyz_xxyy, g_x_0_yyyz_xxyyy, g_x_0_yyyz_xxyyz, g_x_0_yyyz_xxyz, g_x_0_yyyz_xxyzz, g_x_0_yyyz_xxzz, g_x_0_yyyz_xyyy, g_x_0_yyyz_xyyyy, g_x_0_yyyz_xyyyz, g_x_0_yyyz_xyyz, g_x_0_yyyz_xyyzz, g_x_0_yyyz_xyzz, g_x_0_yyyz_xyzzz, g_x_0_yyyz_xzzz, g_x_0_yyyz_yyyy, g_x_0_yyyz_yyyyy, g_x_0_yyyz_yyyyz, g_x_0_yyyz_yyyz, g_x_0_yyyz_yyyzz, g_x_0_yyyz_yyzz, g_x_0_yyyz_yyzzz, g_x_0_yyyz_yzzz, g_x_0_yyyz_yzzzz, g_x_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_xxxx[k] = -g_x_0_yyyz_xxxx[k] * cd_y[k] + g_x_0_yyyz_xxxxy[k];

                g_x_0_yyyyz_xxxy[k] = -g_x_0_yyyz_xxxy[k] * cd_y[k] + g_x_0_yyyz_xxxyy[k];

                g_x_0_yyyyz_xxxz[k] = -g_x_0_yyyz_xxxz[k] * cd_y[k] + g_x_0_yyyz_xxxyz[k];

                g_x_0_yyyyz_xxyy[k] = -g_x_0_yyyz_xxyy[k] * cd_y[k] + g_x_0_yyyz_xxyyy[k];

                g_x_0_yyyyz_xxyz[k] = -g_x_0_yyyz_xxyz[k] * cd_y[k] + g_x_0_yyyz_xxyyz[k];

                g_x_0_yyyyz_xxzz[k] = -g_x_0_yyyz_xxzz[k] * cd_y[k] + g_x_0_yyyz_xxyzz[k];

                g_x_0_yyyyz_xyyy[k] = -g_x_0_yyyz_xyyy[k] * cd_y[k] + g_x_0_yyyz_xyyyy[k];

                g_x_0_yyyyz_xyyz[k] = -g_x_0_yyyz_xyyz[k] * cd_y[k] + g_x_0_yyyz_xyyyz[k];

                g_x_0_yyyyz_xyzz[k] = -g_x_0_yyyz_xyzz[k] * cd_y[k] + g_x_0_yyyz_xyyzz[k];

                g_x_0_yyyyz_xzzz[k] = -g_x_0_yyyz_xzzz[k] * cd_y[k] + g_x_0_yyyz_xyzzz[k];

                g_x_0_yyyyz_yyyy[k] = -g_x_0_yyyz_yyyy[k] * cd_y[k] + g_x_0_yyyz_yyyyy[k];

                g_x_0_yyyyz_yyyz[k] = -g_x_0_yyyz_yyyz[k] * cd_y[k] + g_x_0_yyyz_yyyyz[k];

                g_x_0_yyyyz_yyzz[k] = -g_x_0_yyyz_yyzz[k] * cd_y[k] + g_x_0_yyyz_yyyzz[k];

                g_x_0_yyyyz_yzzz[k] = -g_x_0_yyyz_yzzz[k] * cd_y[k] + g_x_0_yyyz_yyzzz[k];

                g_x_0_yyyyz_zzzz[k] = -g_x_0_yyyz_zzzz[k] * cd_y[k] + g_x_0_yyyz_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_yyyzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_yyyzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_yyyzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_yyyzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_yyyzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_yyyzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_yyyzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_yyyzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_yyyzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_yyyzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_yyyzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_yyyzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_yyyzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_yyyzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 269);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzz_xxxx, g_x_0_yyyzz_xxxy, g_x_0_yyyzz_xxxz, g_x_0_yyyzz_xxyy, g_x_0_yyyzz_xxyz, g_x_0_yyyzz_xxzz, g_x_0_yyyzz_xyyy, g_x_0_yyyzz_xyyz, g_x_0_yyyzz_xyzz, g_x_0_yyyzz_xzzz, g_x_0_yyyzz_yyyy, g_x_0_yyyzz_yyyz, g_x_0_yyyzz_yyzz, g_x_0_yyyzz_yzzz, g_x_0_yyyzz_zzzz, g_x_0_yyzz_xxxx, g_x_0_yyzz_xxxxy, g_x_0_yyzz_xxxy, g_x_0_yyzz_xxxyy, g_x_0_yyzz_xxxyz, g_x_0_yyzz_xxxz, g_x_0_yyzz_xxyy, g_x_0_yyzz_xxyyy, g_x_0_yyzz_xxyyz, g_x_0_yyzz_xxyz, g_x_0_yyzz_xxyzz, g_x_0_yyzz_xxzz, g_x_0_yyzz_xyyy, g_x_0_yyzz_xyyyy, g_x_0_yyzz_xyyyz, g_x_0_yyzz_xyyz, g_x_0_yyzz_xyyzz, g_x_0_yyzz_xyzz, g_x_0_yyzz_xyzzz, g_x_0_yyzz_xzzz, g_x_0_yyzz_yyyy, g_x_0_yyzz_yyyyy, g_x_0_yyzz_yyyyz, g_x_0_yyzz_yyyz, g_x_0_yyzz_yyyzz, g_x_0_yyzz_yyzz, g_x_0_yyzz_yyzzz, g_x_0_yyzz_yzzz, g_x_0_yyzz_yzzzz, g_x_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_xxxx[k] = -g_x_0_yyzz_xxxx[k] * cd_y[k] + g_x_0_yyzz_xxxxy[k];

                g_x_0_yyyzz_xxxy[k] = -g_x_0_yyzz_xxxy[k] * cd_y[k] + g_x_0_yyzz_xxxyy[k];

                g_x_0_yyyzz_xxxz[k] = -g_x_0_yyzz_xxxz[k] * cd_y[k] + g_x_0_yyzz_xxxyz[k];

                g_x_0_yyyzz_xxyy[k] = -g_x_0_yyzz_xxyy[k] * cd_y[k] + g_x_0_yyzz_xxyyy[k];

                g_x_0_yyyzz_xxyz[k] = -g_x_0_yyzz_xxyz[k] * cd_y[k] + g_x_0_yyzz_xxyyz[k];

                g_x_0_yyyzz_xxzz[k] = -g_x_0_yyzz_xxzz[k] * cd_y[k] + g_x_0_yyzz_xxyzz[k];

                g_x_0_yyyzz_xyyy[k] = -g_x_0_yyzz_xyyy[k] * cd_y[k] + g_x_0_yyzz_xyyyy[k];

                g_x_0_yyyzz_xyyz[k] = -g_x_0_yyzz_xyyz[k] * cd_y[k] + g_x_0_yyzz_xyyyz[k];

                g_x_0_yyyzz_xyzz[k] = -g_x_0_yyzz_xyzz[k] * cd_y[k] + g_x_0_yyzz_xyyzz[k];

                g_x_0_yyyzz_xzzz[k] = -g_x_0_yyzz_xzzz[k] * cd_y[k] + g_x_0_yyzz_xyzzz[k];

                g_x_0_yyyzz_yyyy[k] = -g_x_0_yyzz_yyyy[k] * cd_y[k] + g_x_0_yyzz_yyyyy[k];

                g_x_0_yyyzz_yyyz[k] = -g_x_0_yyzz_yyyz[k] * cd_y[k] + g_x_0_yyzz_yyyyz[k];

                g_x_0_yyyzz_yyzz[k] = -g_x_0_yyzz_yyzz[k] * cd_y[k] + g_x_0_yyzz_yyyzz[k];

                g_x_0_yyyzz_yzzz[k] = -g_x_0_yyzz_yzzz[k] * cd_y[k] + g_x_0_yyzz_yyzzz[k];

                g_x_0_yyyzz_zzzz[k] = -g_x_0_yyzz_zzzz[k] * cd_y[k] + g_x_0_yyzz_yzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_yyzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_yyzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_yyzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_yyzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_yyzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_yyzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_yyzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_yyzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_yyzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_yyzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_yyzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_yyzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_yyzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_yyzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 284);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzz_xxxx, g_x_0_yyzzz_xxxy, g_x_0_yyzzz_xxxz, g_x_0_yyzzz_xxyy, g_x_0_yyzzz_xxyz, g_x_0_yyzzz_xxzz, g_x_0_yyzzz_xyyy, g_x_0_yyzzz_xyyz, g_x_0_yyzzz_xyzz, g_x_0_yyzzz_xzzz, g_x_0_yyzzz_yyyy, g_x_0_yyzzz_yyyz, g_x_0_yyzzz_yyzz, g_x_0_yyzzz_yzzz, g_x_0_yyzzz_zzzz, g_x_0_yzzz_xxxx, g_x_0_yzzz_xxxxy, g_x_0_yzzz_xxxy, g_x_0_yzzz_xxxyy, g_x_0_yzzz_xxxyz, g_x_0_yzzz_xxxz, g_x_0_yzzz_xxyy, g_x_0_yzzz_xxyyy, g_x_0_yzzz_xxyyz, g_x_0_yzzz_xxyz, g_x_0_yzzz_xxyzz, g_x_0_yzzz_xxzz, g_x_0_yzzz_xyyy, g_x_0_yzzz_xyyyy, g_x_0_yzzz_xyyyz, g_x_0_yzzz_xyyz, g_x_0_yzzz_xyyzz, g_x_0_yzzz_xyzz, g_x_0_yzzz_xyzzz, g_x_0_yzzz_xzzz, g_x_0_yzzz_yyyy, g_x_0_yzzz_yyyyy, g_x_0_yzzz_yyyyz, g_x_0_yzzz_yyyz, g_x_0_yzzz_yyyzz, g_x_0_yzzz_yyzz, g_x_0_yzzz_yyzzz, g_x_0_yzzz_yzzz, g_x_0_yzzz_yzzzz, g_x_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_xxxx[k] = -g_x_0_yzzz_xxxx[k] * cd_y[k] + g_x_0_yzzz_xxxxy[k];

                g_x_0_yyzzz_xxxy[k] = -g_x_0_yzzz_xxxy[k] * cd_y[k] + g_x_0_yzzz_xxxyy[k];

                g_x_0_yyzzz_xxxz[k] = -g_x_0_yzzz_xxxz[k] * cd_y[k] + g_x_0_yzzz_xxxyz[k];

                g_x_0_yyzzz_xxyy[k] = -g_x_0_yzzz_xxyy[k] * cd_y[k] + g_x_0_yzzz_xxyyy[k];

                g_x_0_yyzzz_xxyz[k] = -g_x_0_yzzz_xxyz[k] * cd_y[k] + g_x_0_yzzz_xxyyz[k];

                g_x_0_yyzzz_xxzz[k] = -g_x_0_yzzz_xxzz[k] * cd_y[k] + g_x_0_yzzz_xxyzz[k];

                g_x_0_yyzzz_xyyy[k] = -g_x_0_yzzz_xyyy[k] * cd_y[k] + g_x_0_yzzz_xyyyy[k];

                g_x_0_yyzzz_xyyz[k] = -g_x_0_yzzz_xyyz[k] * cd_y[k] + g_x_0_yzzz_xyyyz[k];

                g_x_0_yyzzz_xyzz[k] = -g_x_0_yzzz_xyzz[k] * cd_y[k] + g_x_0_yzzz_xyyzz[k];

                g_x_0_yyzzz_xzzz[k] = -g_x_0_yzzz_xzzz[k] * cd_y[k] + g_x_0_yzzz_xyzzz[k];

                g_x_0_yyzzz_yyyy[k] = -g_x_0_yzzz_yyyy[k] * cd_y[k] + g_x_0_yzzz_yyyyy[k];

                g_x_0_yyzzz_yyyz[k] = -g_x_0_yzzz_yyyz[k] * cd_y[k] + g_x_0_yzzz_yyyyz[k];

                g_x_0_yyzzz_yyzz[k] = -g_x_0_yzzz_yyzz[k] * cd_y[k] + g_x_0_yzzz_yyyzz[k];

                g_x_0_yyzzz_yzzz[k] = -g_x_0_yzzz_yzzz[k] * cd_y[k] + g_x_0_yzzz_yyzzz[k];

                g_x_0_yyzzz_zzzz[k] = -g_x_0_yzzz_zzzz[k] * cd_y[k] + g_x_0_yzzz_yzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_yzzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_yzzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_yzzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_yzzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_yzzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_yzzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_yzzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_yzzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_yzzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_yzzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_yzzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_yzzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_yzzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_yzzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 299);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzz_xxxx, g_x_0_yzzzz_xxxy, g_x_0_yzzzz_xxxz, g_x_0_yzzzz_xxyy, g_x_0_yzzzz_xxyz, g_x_0_yzzzz_xxzz, g_x_0_yzzzz_xyyy, g_x_0_yzzzz_xyyz, g_x_0_yzzzz_xyzz, g_x_0_yzzzz_xzzz, g_x_0_yzzzz_yyyy, g_x_0_yzzzz_yyyz, g_x_0_yzzzz_yyzz, g_x_0_yzzzz_yzzz, g_x_0_yzzzz_zzzz, g_x_0_zzzz_xxxx, g_x_0_zzzz_xxxxy, g_x_0_zzzz_xxxy, g_x_0_zzzz_xxxyy, g_x_0_zzzz_xxxyz, g_x_0_zzzz_xxxz, g_x_0_zzzz_xxyy, g_x_0_zzzz_xxyyy, g_x_0_zzzz_xxyyz, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxyzz, g_x_0_zzzz_xxzz, g_x_0_zzzz_xyyy, g_x_0_zzzz_xyyyy, g_x_0_zzzz_xyyyz, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyyzz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xyzzz, g_x_0_zzzz_xzzz, g_x_0_zzzz_yyyy, g_x_0_zzzz_yyyyy, g_x_0_zzzz_yyyyz, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyyzz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yyzzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_yzzzz, g_x_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_xxxx[k] = -g_x_0_zzzz_xxxx[k] * cd_y[k] + g_x_0_zzzz_xxxxy[k];

                g_x_0_yzzzz_xxxy[k] = -g_x_0_zzzz_xxxy[k] * cd_y[k] + g_x_0_zzzz_xxxyy[k];

                g_x_0_yzzzz_xxxz[k] = -g_x_0_zzzz_xxxz[k] * cd_y[k] + g_x_0_zzzz_xxxyz[k];

                g_x_0_yzzzz_xxyy[k] = -g_x_0_zzzz_xxyy[k] * cd_y[k] + g_x_0_zzzz_xxyyy[k];

                g_x_0_yzzzz_xxyz[k] = -g_x_0_zzzz_xxyz[k] * cd_y[k] + g_x_0_zzzz_xxyyz[k];

                g_x_0_yzzzz_xxzz[k] = -g_x_0_zzzz_xxzz[k] * cd_y[k] + g_x_0_zzzz_xxyzz[k];

                g_x_0_yzzzz_xyyy[k] = -g_x_0_zzzz_xyyy[k] * cd_y[k] + g_x_0_zzzz_xyyyy[k];

                g_x_0_yzzzz_xyyz[k] = -g_x_0_zzzz_xyyz[k] * cd_y[k] + g_x_0_zzzz_xyyyz[k];

                g_x_0_yzzzz_xyzz[k] = -g_x_0_zzzz_xyzz[k] * cd_y[k] + g_x_0_zzzz_xyyzz[k];

                g_x_0_yzzzz_xzzz[k] = -g_x_0_zzzz_xzzz[k] * cd_y[k] + g_x_0_zzzz_xyzzz[k];

                g_x_0_yzzzz_yyyy[k] = -g_x_0_zzzz_yyyy[k] * cd_y[k] + g_x_0_zzzz_yyyyy[k];

                g_x_0_yzzzz_yyyz[k] = -g_x_0_zzzz_yyyz[k] * cd_y[k] + g_x_0_zzzz_yyyyz[k];

                g_x_0_yzzzz_yyzz[k] = -g_x_0_zzzz_yyzz[k] * cd_y[k] + g_x_0_zzzz_yyyzz[k];

                g_x_0_yzzzz_yzzz[k] = -g_x_0_zzzz_yzzz[k] * cd_y[k] + g_x_0_zzzz_yyzzz[k];

                g_x_0_yzzzz_zzzz[k] = -g_x_0_zzzz_zzzz[k] * cd_y[k] + g_x_0_zzzz_yzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_zzzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_zzzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_zzzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_zzzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_zzzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_zzzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_zzzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_zzzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_zzzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_zzzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_zzzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_zzzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_zzzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_zzzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_z, g_x_0_zzzz_xxxx, g_x_0_zzzz_xxxxz, g_x_0_zzzz_xxxy, g_x_0_zzzz_xxxyz, g_x_0_zzzz_xxxz, g_x_0_zzzz_xxxzz, g_x_0_zzzz_xxyy, g_x_0_zzzz_xxyyz, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxyzz, g_x_0_zzzz_xxzz, g_x_0_zzzz_xxzzz, g_x_0_zzzz_xyyy, g_x_0_zzzz_xyyyz, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyyzz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xyzzz, g_x_0_zzzz_xzzz, g_x_0_zzzz_xzzzz, g_x_0_zzzz_yyyy, g_x_0_zzzz_yyyyz, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyyzz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yyzzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_yzzzz, g_x_0_zzzz_zzzz, g_x_0_zzzz_zzzzz, g_x_0_zzzzz_xxxx, g_x_0_zzzzz_xxxy, g_x_0_zzzzz_xxxz, g_x_0_zzzzz_xxyy, g_x_0_zzzzz_xxyz, g_x_0_zzzzz_xxzz, g_x_0_zzzzz_xyyy, g_x_0_zzzzz_xyyz, g_x_0_zzzzz_xyzz, g_x_0_zzzzz_xzzz, g_x_0_zzzzz_yyyy, g_x_0_zzzzz_yyyz, g_x_0_zzzzz_yyzz, g_x_0_zzzzz_yzzz, g_x_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_xxxx[k] = -g_x_0_zzzz_xxxx[k] * cd_z[k] + g_x_0_zzzz_xxxxz[k];

                g_x_0_zzzzz_xxxy[k] = -g_x_0_zzzz_xxxy[k] * cd_z[k] + g_x_0_zzzz_xxxyz[k];

                g_x_0_zzzzz_xxxz[k] = -g_x_0_zzzz_xxxz[k] * cd_z[k] + g_x_0_zzzz_xxxzz[k];

                g_x_0_zzzzz_xxyy[k] = -g_x_0_zzzz_xxyy[k] * cd_z[k] + g_x_0_zzzz_xxyyz[k];

                g_x_0_zzzzz_xxyz[k] = -g_x_0_zzzz_xxyz[k] * cd_z[k] + g_x_0_zzzz_xxyzz[k];

                g_x_0_zzzzz_xxzz[k] = -g_x_0_zzzz_xxzz[k] * cd_z[k] + g_x_0_zzzz_xxzzz[k];

                g_x_0_zzzzz_xyyy[k] = -g_x_0_zzzz_xyyy[k] * cd_z[k] + g_x_0_zzzz_xyyyz[k];

                g_x_0_zzzzz_xyyz[k] = -g_x_0_zzzz_xyyz[k] * cd_z[k] + g_x_0_zzzz_xyyzz[k];

                g_x_0_zzzzz_xyzz[k] = -g_x_0_zzzz_xyzz[k] * cd_z[k] + g_x_0_zzzz_xyzzz[k];

                g_x_0_zzzzz_xzzz[k] = -g_x_0_zzzz_xzzz[k] * cd_z[k] + g_x_0_zzzz_xzzzz[k];

                g_x_0_zzzzz_yyyy[k] = -g_x_0_zzzz_yyyy[k] * cd_z[k] + g_x_0_zzzz_yyyyz[k];

                g_x_0_zzzzz_yyyz[k] = -g_x_0_zzzz_yyyz[k] * cd_z[k] + g_x_0_zzzz_yyyzz[k];

                g_x_0_zzzzz_yyzz[k] = -g_x_0_zzzz_yyzz[k] * cd_z[k] + g_x_0_zzzz_yyzzz[k];

                g_x_0_zzzzz_yzzz[k] = -g_x_0_zzzz_yzzz[k] * cd_z[k] + g_x_0_zzzz_yzzzz[k];

                g_x_0_zzzzz_zzzz[k] = -g_x_0_zzzz_zzzz[k] * cd_z[k] + g_x_0_zzzz_zzzzz[k];
            }
            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxx_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 5);

            auto g_y_0_xxxxx_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 6);

            auto g_y_0_xxxxx_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 7);

            auto g_y_0_xxxxx_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 8);

            auto g_y_0_xxxxx_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 9);

            auto g_y_0_xxxxx_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 10);

            auto g_y_0_xxxxx_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 11);

            auto g_y_0_xxxxx_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 12);

            auto g_y_0_xxxxx_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 13);

            auto g_y_0_xxxxx_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_y_0_xxxx_xxxx, g_y_0_xxxx_xxxxx, g_y_0_xxxx_xxxxy, g_y_0_xxxx_xxxxz, g_y_0_xxxx_xxxy, g_y_0_xxxx_xxxyy, g_y_0_xxxx_xxxyz, g_y_0_xxxx_xxxz, g_y_0_xxxx_xxxzz, g_y_0_xxxx_xxyy, g_y_0_xxxx_xxyyy, g_y_0_xxxx_xxyyz, g_y_0_xxxx_xxyz, g_y_0_xxxx_xxyzz, g_y_0_xxxx_xxzz, g_y_0_xxxx_xxzzz, g_y_0_xxxx_xyyy, g_y_0_xxxx_xyyyy, g_y_0_xxxx_xyyyz, g_y_0_xxxx_xyyz, g_y_0_xxxx_xyyzz, g_y_0_xxxx_xyzz, g_y_0_xxxx_xyzzz, g_y_0_xxxx_xzzz, g_y_0_xxxx_xzzzz, g_y_0_xxxx_yyyy, g_y_0_xxxx_yyyz, g_y_0_xxxx_yyzz, g_y_0_xxxx_yzzz, g_y_0_xxxx_zzzz, g_y_0_xxxxx_xxxx, g_y_0_xxxxx_xxxy, g_y_0_xxxxx_xxxz, g_y_0_xxxxx_xxyy, g_y_0_xxxxx_xxyz, g_y_0_xxxxx_xxzz, g_y_0_xxxxx_xyyy, g_y_0_xxxxx_xyyz, g_y_0_xxxxx_xyzz, g_y_0_xxxxx_xzzz, g_y_0_xxxxx_yyyy, g_y_0_xxxxx_yyyz, g_y_0_xxxxx_yyzz, g_y_0_xxxxx_yzzz, g_y_0_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_xxxx[k] = -g_y_0_xxxx_xxxx[k] * cd_x[k] + g_y_0_xxxx_xxxxx[k];

                g_y_0_xxxxx_xxxy[k] = -g_y_0_xxxx_xxxy[k] * cd_x[k] + g_y_0_xxxx_xxxxy[k];

                g_y_0_xxxxx_xxxz[k] = -g_y_0_xxxx_xxxz[k] * cd_x[k] + g_y_0_xxxx_xxxxz[k];

                g_y_0_xxxxx_xxyy[k] = -g_y_0_xxxx_xxyy[k] * cd_x[k] + g_y_0_xxxx_xxxyy[k];

                g_y_0_xxxxx_xxyz[k] = -g_y_0_xxxx_xxyz[k] * cd_x[k] + g_y_0_xxxx_xxxyz[k];

                g_y_0_xxxxx_xxzz[k] = -g_y_0_xxxx_xxzz[k] * cd_x[k] + g_y_0_xxxx_xxxzz[k];

                g_y_0_xxxxx_xyyy[k] = -g_y_0_xxxx_xyyy[k] * cd_x[k] + g_y_0_xxxx_xxyyy[k];

                g_y_0_xxxxx_xyyz[k] = -g_y_0_xxxx_xyyz[k] * cd_x[k] + g_y_0_xxxx_xxyyz[k];

                g_y_0_xxxxx_xyzz[k] = -g_y_0_xxxx_xyzz[k] * cd_x[k] + g_y_0_xxxx_xxyzz[k];

                g_y_0_xxxxx_xzzz[k] = -g_y_0_xxxx_xzzz[k] * cd_x[k] + g_y_0_xxxx_xxzzz[k];

                g_y_0_xxxxx_yyyy[k] = -g_y_0_xxxx_yyyy[k] * cd_x[k] + g_y_0_xxxx_xyyyy[k];

                g_y_0_xxxxx_yyyz[k] = -g_y_0_xxxx_yyyz[k] * cd_x[k] + g_y_0_xxxx_xyyyz[k];

                g_y_0_xxxxx_yyzz[k] = -g_y_0_xxxx_yyzz[k] * cd_x[k] + g_y_0_xxxx_xyyzz[k];

                g_y_0_xxxxx_yzzz[k] = -g_y_0_xxxx_yzzz[k] * cd_x[k] + g_y_0_xxxx_xyzzz[k];

                g_y_0_xxxxx_zzzz[k] = -g_y_0_xxxx_zzzz[k] * cd_x[k] + g_y_0_xxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxy_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 15);

            auto g_y_0_xxxxy_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 16);

            auto g_y_0_xxxxy_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 17);

            auto g_y_0_xxxxy_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 18);

            auto g_y_0_xxxxy_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 19);

            auto g_y_0_xxxxy_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 20);

            auto g_y_0_xxxxy_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 21);

            auto g_y_0_xxxxy_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 22);

            auto g_y_0_xxxxy_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 23);

            auto g_y_0_xxxxy_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 24);

            auto g_y_0_xxxxy_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 25);

            auto g_y_0_xxxxy_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 26);

            auto g_y_0_xxxxy_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 27);

            auto g_y_0_xxxxy_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 28);

            auto g_y_0_xxxxy_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxy_xxxx, g_y_0_xxxxy_xxxy, g_y_0_xxxxy_xxxz, g_y_0_xxxxy_xxyy, g_y_0_xxxxy_xxyz, g_y_0_xxxxy_xxzz, g_y_0_xxxxy_xyyy, g_y_0_xxxxy_xyyz, g_y_0_xxxxy_xyzz, g_y_0_xxxxy_xzzz, g_y_0_xxxxy_yyyy, g_y_0_xxxxy_yyyz, g_y_0_xxxxy_yyzz, g_y_0_xxxxy_yzzz, g_y_0_xxxxy_zzzz, g_y_0_xxxy_xxxx, g_y_0_xxxy_xxxxx, g_y_0_xxxy_xxxxy, g_y_0_xxxy_xxxxz, g_y_0_xxxy_xxxy, g_y_0_xxxy_xxxyy, g_y_0_xxxy_xxxyz, g_y_0_xxxy_xxxz, g_y_0_xxxy_xxxzz, g_y_0_xxxy_xxyy, g_y_0_xxxy_xxyyy, g_y_0_xxxy_xxyyz, g_y_0_xxxy_xxyz, g_y_0_xxxy_xxyzz, g_y_0_xxxy_xxzz, g_y_0_xxxy_xxzzz, g_y_0_xxxy_xyyy, g_y_0_xxxy_xyyyy, g_y_0_xxxy_xyyyz, g_y_0_xxxy_xyyz, g_y_0_xxxy_xyyzz, g_y_0_xxxy_xyzz, g_y_0_xxxy_xyzzz, g_y_0_xxxy_xzzz, g_y_0_xxxy_xzzzz, g_y_0_xxxy_yyyy, g_y_0_xxxy_yyyz, g_y_0_xxxy_yyzz, g_y_0_xxxy_yzzz, g_y_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_xxxx[k] = -g_y_0_xxxy_xxxx[k] * cd_x[k] + g_y_0_xxxy_xxxxx[k];

                g_y_0_xxxxy_xxxy[k] = -g_y_0_xxxy_xxxy[k] * cd_x[k] + g_y_0_xxxy_xxxxy[k];

                g_y_0_xxxxy_xxxz[k] = -g_y_0_xxxy_xxxz[k] * cd_x[k] + g_y_0_xxxy_xxxxz[k];

                g_y_0_xxxxy_xxyy[k] = -g_y_0_xxxy_xxyy[k] * cd_x[k] + g_y_0_xxxy_xxxyy[k];

                g_y_0_xxxxy_xxyz[k] = -g_y_0_xxxy_xxyz[k] * cd_x[k] + g_y_0_xxxy_xxxyz[k];

                g_y_0_xxxxy_xxzz[k] = -g_y_0_xxxy_xxzz[k] * cd_x[k] + g_y_0_xxxy_xxxzz[k];

                g_y_0_xxxxy_xyyy[k] = -g_y_0_xxxy_xyyy[k] * cd_x[k] + g_y_0_xxxy_xxyyy[k];

                g_y_0_xxxxy_xyyz[k] = -g_y_0_xxxy_xyyz[k] * cd_x[k] + g_y_0_xxxy_xxyyz[k];

                g_y_0_xxxxy_xyzz[k] = -g_y_0_xxxy_xyzz[k] * cd_x[k] + g_y_0_xxxy_xxyzz[k];

                g_y_0_xxxxy_xzzz[k] = -g_y_0_xxxy_xzzz[k] * cd_x[k] + g_y_0_xxxy_xxzzz[k];

                g_y_0_xxxxy_yyyy[k] = -g_y_0_xxxy_yyyy[k] * cd_x[k] + g_y_0_xxxy_xyyyy[k];

                g_y_0_xxxxy_yyyz[k] = -g_y_0_xxxy_yyyz[k] * cd_x[k] + g_y_0_xxxy_xyyyz[k];

                g_y_0_xxxxy_yyzz[k] = -g_y_0_xxxy_yyzz[k] * cd_x[k] + g_y_0_xxxy_xyyzz[k];

                g_y_0_xxxxy_yzzz[k] = -g_y_0_xxxy_yzzz[k] * cd_x[k] + g_y_0_xxxy_xyzzz[k];

                g_y_0_xxxxy_zzzz[k] = -g_y_0_xxxy_zzzz[k] * cd_x[k] + g_y_0_xxxy_xzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 30);

            auto g_y_0_xxxxz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 31);

            auto g_y_0_xxxxz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 32);

            auto g_y_0_xxxxz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 33);

            auto g_y_0_xxxxz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 34);

            auto g_y_0_xxxxz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 35);

            auto g_y_0_xxxxz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 36);

            auto g_y_0_xxxxz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 37);

            auto g_y_0_xxxxz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 38);

            auto g_y_0_xxxxz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 39);

            auto g_y_0_xxxxz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 40);

            auto g_y_0_xxxxz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 41);

            auto g_y_0_xxxxz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 42);

            auto g_y_0_xxxxz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 43);

            auto g_y_0_xxxxz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxz_xxxx, g_y_0_xxxxz_xxxy, g_y_0_xxxxz_xxxz, g_y_0_xxxxz_xxyy, g_y_0_xxxxz_xxyz, g_y_0_xxxxz_xxzz, g_y_0_xxxxz_xyyy, g_y_0_xxxxz_xyyz, g_y_0_xxxxz_xyzz, g_y_0_xxxxz_xzzz, g_y_0_xxxxz_yyyy, g_y_0_xxxxz_yyyz, g_y_0_xxxxz_yyzz, g_y_0_xxxxz_yzzz, g_y_0_xxxxz_zzzz, g_y_0_xxxz_xxxx, g_y_0_xxxz_xxxxx, g_y_0_xxxz_xxxxy, g_y_0_xxxz_xxxxz, g_y_0_xxxz_xxxy, g_y_0_xxxz_xxxyy, g_y_0_xxxz_xxxyz, g_y_0_xxxz_xxxz, g_y_0_xxxz_xxxzz, g_y_0_xxxz_xxyy, g_y_0_xxxz_xxyyy, g_y_0_xxxz_xxyyz, g_y_0_xxxz_xxyz, g_y_0_xxxz_xxyzz, g_y_0_xxxz_xxzz, g_y_0_xxxz_xxzzz, g_y_0_xxxz_xyyy, g_y_0_xxxz_xyyyy, g_y_0_xxxz_xyyyz, g_y_0_xxxz_xyyz, g_y_0_xxxz_xyyzz, g_y_0_xxxz_xyzz, g_y_0_xxxz_xyzzz, g_y_0_xxxz_xzzz, g_y_0_xxxz_xzzzz, g_y_0_xxxz_yyyy, g_y_0_xxxz_yyyz, g_y_0_xxxz_yyzz, g_y_0_xxxz_yzzz, g_y_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_xxxx[k] = -g_y_0_xxxz_xxxx[k] * cd_x[k] + g_y_0_xxxz_xxxxx[k];

                g_y_0_xxxxz_xxxy[k] = -g_y_0_xxxz_xxxy[k] * cd_x[k] + g_y_0_xxxz_xxxxy[k];

                g_y_0_xxxxz_xxxz[k] = -g_y_0_xxxz_xxxz[k] * cd_x[k] + g_y_0_xxxz_xxxxz[k];

                g_y_0_xxxxz_xxyy[k] = -g_y_0_xxxz_xxyy[k] * cd_x[k] + g_y_0_xxxz_xxxyy[k];

                g_y_0_xxxxz_xxyz[k] = -g_y_0_xxxz_xxyz[k] * cd_x[k] + g_y_0_xxxz_xxxyz[k];

                g_y_0_xxxxz_xxzz[k] = -g_y_0_xxxz_xxzz[k] * cd_x[k] + g_y_0_xxxz_xxxzz[k];

                g_y_0_xxxxz_xyyy[k] = -g_y_0_xxxz_xyyy[k] * cd_x[k] + g_y_0_xxxz_xxyyy[k];

                g_y_0_xxxxz_xyyz[k] = -g_y_0_xxxz_xyyz[k] * cd_x[k] + g_y_0_xxxz_xxyyz[k];

                g_y_0_xxxxz_xyzz[k] = -g_y_0_xxxz_xyzz[k] * cd_x[k] + g_y_0_xxxz_xxyzz[k];

                g_y_0_xxxxz_xzzz[k] = -g_y_0_xxxz_xzzz[k] * cd_x[k] + g_y_0_xxxz_xxzzz[k];

                g_y_0_xxxxz_yyyy[k] = -g_y_0_xxxz_yyyy[k] * cd_x[k] + g_y_0_xxxz_xyyyy[k];

                g_y_0_xxxxz_yyyz[k] = -g_y_0_xxxz_yyyz[k] * cd_x[k] + g_y_0_xxxz_xyyyz[k];

                g_y_0_xxxxz_yyzz[k] = -g_y_0_xxxz_yyzz[k] * cd_x[k] + g_y_0_xxxz_xyyzz[k];

                g_y_0_xxxxz_yzzz[k] = -g_y_0_xxxz_yzzz[k] * cd_x[k] + g_y_0_xxxz_xyzzz[k];

                g_y_0_xxxxz_zzzz[k] = -g_y_0_xxxz_zzzz[k] * cd_x[k] + g_y_0_xxxz_xzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyy_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 45);

            auto g_y_0_xxxyy_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 46);

            auto g_y_0_xxxyy_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 47);

            auto g_y_0_xxxyy_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 48);

            auto g_y_0_xxxyy_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 49);

            auto g_y_0_xxxyy_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 50);

            auto g_y_0_xxxyy_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 51);

            auto g_y_0_xxxyy_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 52);

            auto g_y_0_xxxyy_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 53);

            auto g_y_0_xxxyy_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 54);

            auto g_y_0_xxxyy_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 55);

            auto g_y_0_xxxyy_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 56);

            auto g_y_0_xxxyy_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 57);

            auto g_y_0_xxxyy_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 58);

            auto g_y_0_xxxyy_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyy_xxxx, g_y_0_xxxyy_xxxy, g_y_0_xxxyy_xxxz, g_y_0_xxxyy_xxyy, g_y_0_xxxyy_xxyz, g_y_0_xxxyy_xxzz, g_y_0_xxxyy_xyyy, g_y_0_xxxyy_xyyz, g_y_0_xxxyy_xyzz, g_y_0_xxxyy_xzzz, g_y_0_xxxyy_yyyy, g_y_0_xxxyy_yyyz, g_y_0_xxxyy_yyzz, g_y_0_xxxyy_yzzz, g_y_0_xxxyy_zzzz, g_y_0_xxyy_xxxx, g_y_0_xxyy_xxxxx, g_y_0_xxyy_xxxxy, g_y_0_xxyy_xxxxz, g_y_0_xxyy_xxxy, g_y_0_xxyy_xxxyy, g_y_0_xxyy_xxxyz, g_y_0_xxyy_xxxz, g_y_0_xxyy_xxxzz, g_y_0_xxyy_xxyy, g_y_0_xxyy_xxyyy, g_y_0_xxyy_xxyyz, g_y_0_xxyy_xxyz, g_y_0_xxyy_xxyzz, g_y_0_xxyy_xxzz, g_y_0_xxyy_xxzzz, g_y_0_xxyy_xyyy, g_y_0_xxyy_xyyyy, g_y_0_xxyy_xyyyz, g_y_0_xxyy_xyyz, g_y_0_xxyy_xyyzz, g_y_0_xxyy_xyzz, g_y_0_xxyy_xyzzz, g_y_0_xxyy_xzzz, g_y_0_xxyy_xzzzz, g_y_0_xxyy_yyyy, g_y_0_xxyy_yyyz, g_y_0_xxyy_yyzz, g_y_0_xxyy_yzzz, g_y_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_xxxx[k] = -g_y_0_xxyy_xxxx[k] * cd_x[k] + g_y_0_xxyy_xxxxx[k];

                g_y_0_xxxyy_xxxy[k] = -g_y_0_xxyy_xxxy[k] * cd_x[k] + g_y_0_xxyy_xxxxy[k];

                g_y_0_xxxyy_xxxz[k] = -g_y_0_xxyy_xxxz[k] * cd_x[k] + g_y_0_xxyy_xxxxz[k];

                g_y_0_xxxyy_xxyy[k] = -g_y_0_xxyy_xxyy[k] * cd_x[k] + g_y_0_xxyy_xxxyy[k];

                g_y_0_xxxyy_xxyz[k] = -g_y_0_xxyy_xxyz[k] * cd_x[k] + g_y_0_xxyy_xxxyz[k];

                g_y_0_xxxyy_xxzz[k] = -g_y_0_xxyy_xxzz[k] * cd_x[k] + g_y_0_xxyy_xxxzz[k];

                g_y_0_xxxyy_xyyy[k] = -g_y_0_xxyy_xyyy[k] * cd_x[k] + g_y_0_xxyy_xxyyy[k];

                g_y_0_xxxyy_xyyz[k] = -g_y_0_xxyy_xyyz[k] * cd_x[k] + g_y_0_xxyy_xxyyz[k];

                g_y_0_xxxyy_xyzz[k] = -g_y_0_xxyy_xyzz[k] * cd_x[k] + g_y_0_xxyy_xxyzz[k];

                g_y_0_xxxyy_xzzz[k] = -g_y_0_xxyy_xzzz[k] * cd_x[k] + g_y_0_xxyy_xxzzz[k];

                g_y_0_xxxyy_yyyy[k] = -g_y_0_xxyy_yyyy[k] * cd_x[k] + g_y_0_xxyy_xyyyy[k];

                g_y_0_xxxyy_yyyz[k] = -g_y_0_xxyy_yyyz[k] * cd_x[k] + g_y_0_xxyy_xyyyz[k];

                g_y_0_xxxyy_yyzz[k] = -g_y_0_xxyy_yyzz[k] * cd_x[k] + g_y_0_xxyy_xyyzz[k];

                g_y_0_xxxyy_yzzz[k] = -g_y_0_xxyy_yzzz[k] * cd_x[k] + g_y_0_xxyy_xyzzz[k];

                g_y_0_xxxyy_zzzz[k] = -g_y_0_xxyy_zzzz[k] * cd_x[k] + g_y_0_xxyy_xzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 60);

            auto g_y_0_xxxyz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 61);

            auto g_y_0_xxxyz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 62);

            auto g_y_0_xxxyz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 63);

            auto g_y_0_xxxyz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 64);

            auto g_y_0_xxxyz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 65);

            auto g_y_0_xxxyz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 66);

            auto g_y_0_xxxyz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 67);

            auto g_y_0_xxxyz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 68);

            auto g_y_0_xxxyz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 69);

            auto g_y_0_xxxyz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 70);

            auto g_y_0_xxxyz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 71);

            auto g_y_0_xxxyz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 72);

            auto g_y_0_xxxyz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 73);

            auto g_y_0_xxxyz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 74);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyz_xxxx, g_y_0_xxxyz_xxxy, g_y_0_xxxyz_xxxz, g_y_0_xxxyz_xxyy, g_y_0_xxxyz_xxyz, g_y_0_xxxyz_xxzz, g_y_0_xxxyz_xyyy, g_y_0_xxxyz_xyyz, g_y_0_xxxyz_xyzz, g_y_0_xxxyz_xzzz, g_y_0_xxxyz_yyyy, g_y_0_xxxyz_yyyz, g_y_0_xxxyz_yyzz, g_y_0_xxxyz_yzzz, g_y_0_xxxyz_zzzz, g_y_0_xxyz_xxxx, g_y_0_xxyz_xxxxx, g_y_0_xxyz_xxxxy, g_y_0_xxyz_xxxxz, g_y_0_xxyz_xxxy, g_y_0_xxyz_xxxyy, g_y_0_xxyz_xxxyz, g_y_0_xxyz_xxxz, g_y_0_xxyz_xxxzz, g_y_0_xxyz_xxyy, g_y_0_xxyz_xxyyy, g_y_0_xxyz_xxyyz, g_y_0_xxyz_xxyz, g_y_0_xxyz_xxyzz, g_y_0_xxyz_xxzz, g_y_0_xxyz_xxzzz, g_y_0_xxyz_xyyy, g_y_0_xxyz_xyyyy, g_y_0_xxyz_xyyyz, g_y_0_xxyz_xyyz, g_y_0_xxyz_xyyzz, g_y_0_xxyz_xyzz, g_y_0_xxyz_xyzzz, g_y_0_xxyz_xzzz, g_y_0_xxyz_xzzzz, g_y_0_xxyz_yyyy, g_y_0_xxyz_yyyz, g_y_0_xxyz_yyzz, g_y_0_xxyz_yzzz, g_y_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_xxxx[k] = -g_y_0_xxyz_xxxx[k] * cd_x[k] + g_y_0_xxyz_xxxxx[k];

                g_y_0_xxxyz_xxxy[k] = -g_y_0_xxyz_xxxy[k] * cd_x[k] + g_y_0_xxyz_xxxxy[k];

                g_y_0_xxxyz_xxxz[k] = -g_y_0_xxyz_xxxz[k] * cd_x[k] + g_y_0_xxyz_xxxxz[k];

                g_y_0_xxxyz_xxyy[k] = -g_y_0_xxyz_xxyy[k] * cd_x[k] + g_y_0_xxyz_xxxyy[k];

                g_y_0_xxxyz_xxyz[k] = -g_y_0_xxyz_xxyz[k] * cd_x[k] + g_y_0_xxyz_xxxyz[k];

                g_y_0_xxxyz_xxzz[k] = -g_y_0_xxyz_xxzz[k] * cd_x[k] + g_y_0_xxyz_xxxzz[k];

                g_y_0_xxxyz_xyyy[k] = -g_y_0_xxyz_xyyy[k] * cd_x[k] + g_y_0_xxyz_xxyyy[k];

                g_y_0_xxxyz_xyyz[k] = -g_y_0_xxyz_xyyz[k] * cd_x[k] + g_y_0_xxyz_xxyyz[k];

                g_y_0_xxxyz_xyzz[k] = -g_y_0_xxyz_xyzz[k] * cd_x[k] + g_y_0_xxyz_xxyzz[k];

                g_y_0_xxxyz_xzzz[k] = -g_y_0_xxyz_xzzz[k] * cd_x[k] + g_y_0_xxyz_xxzzz[k];

                g_y_0_xxxyz_yyyy[k] = -g_y_0_xxyz_yyyy[k] * cd_x[k] + g_y_0_xxyz_xyyyy[k];

                g_y_0_xxxyz_yyyz[k] = -g_y_0_xxyz_yyyz[k] * cd_x[k] + g_y_0_xxyz_xyyyz[k];

                g_y_0_xxxyz_yyzz[k] = -g_y_0_xxyz_yyzz[k] * cd_x[k] + g_y_0_xxyz_xyyzz[k];

                g_y_0_xxxyz_yzzz[k] = -g_y_0_xxyz_yzzz[k] * cd_x[k] + g_y_0_xxyz_xyzzz[k];

                g_y_0_xxxyz_zzzz[k] = -g_y_0_xxyz_zzzz[k] * cd_x[k] + g_y_0_xxyz_xzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 75);

            auto g_y_0_xxxzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 76);

            auto g_y_0_xxxzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 77);

            auto g_y_0_xxxzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 78);

            auto g_y_0_xxxzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 79);

            auto g_y_0_xxxzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 80);

            auto g_y_0_xxxzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 81);

            auto g_y_0_xxxzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 82);

            auto g_y_0_xxxzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 83);

            auto g_y_0_xxxzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 84);

            auto g_y_0_xxxzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 85);

            auto g_y_0_xxxzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 86);

            auto g_y_0_xxxzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 87);

            auto g_y_0_xxxzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 88);

            auto g_y_0_xxxzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzz_xxxx, g_y_0_xxxzz_xxxy, g_y_0_xxxzz_xxxz, g_y_0_xxxzz_xxyy, g_y_0_xxxzz_xxyz, g_y_0_xxxzz_xxzz, g_y_0_xxxzz_xyyy, g_y_0_xxxzz_xyyz, g_y_0_xxxzz_xyzz, g_y_0_xxxzz_xzzz, g_y_0_xxxzz_yyyy, g_y_0_xxxzz_yyyz, g_y_0_xxxzz_yyzz, g_y_0_xxxzz_yzzz, g_y_0_xxxzz_zzzz, g_y_0_xxzz_xxxx, g_y_0_xxzz_xxxxx, g_y_0_xxzz_xxxxy, g_y_0_xxzz_xxxxz, g_y_0_xxzz_xxxy, g_y_0_xxzz_xxxyy, g_y_0_xxzz_xxxyz, g_y_0_xxzz_xxxz, g_y_0_xxzz_xxxzz, g_y_0_xxzz_xxyy, g_y_0_xxzz_xxyyy, g_y_0_xxzz_xxyyz, g_y_0_xxzz_xxyz, g_y_0_xxzz_xxyzz, g_y_0_xxzz_xxzz, g_y_0_xxzz_xxzzz, g_y_0_xxzz_xyyy, g_y_0_xxzz_xyyyy, g_y_0_xxzz_xyyyz, g_y_0_xxzz_xyyz, g_y_0_xxzz_xyyzz, g_y_0_xxzz_xyzz, g_y_0_xxzz_xyzzz, g_y_0_xxzz_xzzz, g_y_0_xxzz_xzzzz, g_y_0_xxzz_yyyy, g_y_0_xxzz_yyyz, g_y_0_xxzz_yyzz, g_y_0_xxzz_yzzz, g_y_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_xxxx[k] = -g_y_0_xxzz_xxxx[k] * cd_x[k] + g_y_0_xxzz_xxxxx[k];

                g_y_0_xxxzz_xxxy[k] = -g_y_0_xxzz_xxxy[k] * cd_x[k] + g_y_0_xxzz_xxxxy[k];

                g_y_0_xxxzz_xxxz[k] = -g_y_0_xxzz_xxxz[k] * cd_x[k] + g_y_0_xxzz_xxxxz[k];

                g_y_0_xxxzz_xxyy[k] = -g_y_0_xxzz_xxyy[k] * cd_x[k] + g_y_0_xxzz_xxxyy[k];

                g_y_0_xxxzz_xxyz[k] = -g_y_0_xxzz_xxyz[k] * cd_x[k] + g_y_0_xxzz_xxxyz[k];

                g_y_0_xxxzz_xxzz[k] = -g_y_0_xxzz_xxzz[k] * cd_x[k] + g_y_0_xxzz_xxxzz[k];

                g_y_0_xxxzz_xyyy[k] = -g_y_0_xxzz_xyyy[k] * cd_x[k] + g_y_0_xxzz_xxyyy[k];

                g_y_0_xxxzz_xyyz[k] = -g_y_0_xxzz_xyyz[k] * cd_x[k] + g_y_0_xxzz_xxyyz[k];

                g_y_0_xxxzz_xyzz[k] = -g_y_0_xxzz_xyzz[k] * cd_x[k] + g_y_0_xxzz_xxyzz[k];

                g_y_0_xxxzz_xzzz[k] = -g_y_0_xxzz_xzzz[k] * cd_x[k] + g_y_0_xxzz_xxzzz[k];

                g_y_0_xxxzz_yyyy[k] = -g_y_0_xxzz_yyyy[k] * cd_x[k] + g_y_0_xxzz_xyyyy[k];

                g_y_0_xxxzz_yyyz[k] = -g_y_0_xxzz_yyyz[k] * cd_x[k] + g_y_0_xxzz_xyyyz[k];

                g_y_0_xxxzz_yyzz[k] = -g_y_0_xxzz_yyzz[k] * cd_x[k] + g_y_0_xxzz_xyyzz[k];

                g_y_0_xxxzz_yzzz[k] = -g_y_0_xxzz_yzzz[k] * cd_x[k] + g_y_0_xxzz_xyzzz[k];

                g_y_0_xxxzz_zzzz[k] = -g_y_0_xxzz_zzzz[k] * cd_x[k] + g_y_0_xxzz_xzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyy_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 90);

            auto g_y_0_xxyyy_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 91);

            auto g_y_0_xxyyy_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 92);

            auto g_y_0_xxyyy_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 93);

            auto g_y_0_xxyyy_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 94);

            auto g_y_0_xxyyy_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 95);

            auto g_y_0_xxyyy_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 96);

            auto g_y_0_xxyyy_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 97);

            auto g_y_0_xxyyy_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 98);

            auto g_y_0_xxyyy_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 99);

            auto g_y_0_xxyyy_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 100);

            auto g_y_0_xxyyy_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 101);

            auto g_y_0_xxyyy_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 102);

            auto g_y_0_xxyyy_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 103);

            auto g_y_0_xxyyy_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyy_xxxx, g_y_0_xxyyy_xxxy, g_y_0_xxyyy_xxxz, g_y_0_xxyyy_xxyy, g_y_0_xxyyy_xxyz, g_y_0_xxyyy_xxzz, g_y_0_xxyyy_xyyy, g_y_0_xxyyy_xyyz, g_y_0_xxyyy_xyzz, g_y_0_xxyyy_xzzz, g_y_0_xxyyy_yyyy, g_y_0_xxyyy_yyyz, g_y_0_xxyyy_yyzz, g_y_0_xxyyy_yzzz, g_y_0_xxyyy_zzzz, g_y_0_xyyy_xxxx, g_y_0_xyyy_xxxxx, g_y_0_xyyy_xxxxy, g_y_0_xyyy_xxxxz, g_y_0_xyyy_xxxy, g_y_0_xyyy_xxxyy, g_y_0_xyyy_xxxyz, g_y_0_xyyy_xxxz, g_y_0_xyyy_xxxzz, g_y_0_xyyy_xxyy, g_y_0_xyyy_xxyyy, g_y_0_xyyy_xxyyz, g_y_0_xyyy_xxyz, g_y_0_xyyy_xxyzz, g_y_0_xyyy_xxzz, g_y_0_xyyy_xxzzz, g_y_0_xyyy_xyyy, g_y_0_xyyy_xyyyy, g_y_0_xyyy_xyyyz, g_y_0_xyyy_xyyz, g_y_0_xyyy_xyyzz, g_y_0_xyyy_xyzz, g_y_0_xyyy_xyzzz, g_y_0_xyyy_xzzz, g_y_0_xyyy_xzzzz, g_y_0_xyyy_yyyy, g_y_0_xyyy_yyyz, g_y_0_xyyy_yyzz, g_y_0_xyyy_yzzz, g_y_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_xxxx[k] = -g_y_0_xyyy_xxxx[k] * cd_x[k] + g_y_0_xyyy_xxxxx[k];

                g_y_0_xxyyy_xxxy[k] = -g_y_0_xyyy_xxxy[k] * cd_x[k] + g_y_0_xyyy_xxxxy[k];

                g_y_0_xxyyy_xxxz[k] = -g_y_0_xyyy_xxxz[k] * cd_x[k] + g_y_0_xyyy_xxxxz[k];

                g_y_0_xxyyy_xxyy[k] = -g_y_0_xyyy_xxyy[k] * cd_x[k] + g_y_0_xyyy_xxxyy[k];

                g_y_0_xxyyy_xxyz[k] = -g_y_0_xyyy_xxyz[k] * cd_x[k] + g_y_0_xyyy_xxxyz[k];

                g_y_0_xxyyy_xxzz[k] = -g_y_0_xyyy_xxzz[k] * cd_x[k] + g_y_0_xyyy_xxxzz[k];

                g_y_0_xxyyy_xyyy[k] = -g_y_0_xyyy_xyyy[k] * cd_x[k] + g_y_0_xyyy_xxyyy[k];

                g_y_0_xxyyy_xyyz[k] = -g_y_0_xyyy_xyyz[k] * cd_x[k] + g_y_0_xyyy_xxyyz[k];

                g_y_0_xxyyy_xyzz[k] = -g_y_0_xyyy_xyzz[k] * cd_x[k] + g_y_0_xyyy_xxyzz[k];

                g_y_0_xxyyy_xzzz[k] = -g_y_0_xyyy_xzzz[k] * cd_x[k] + g_y_0_xyyy_xxzzz[k];

                g_y_0_xxyyy_yyyy[k] = -g_y_0_xyyy_yyyy[k] * cd_x[k] + g_y_0_xyyy_xyyyy[k];

                g_y_0_xxyyy_yyyz[k] = -g_y_0_xyyy_yyyz[k] * cd_x[k] + g_y_0_xyyy_xyyyz[k];

                g_y_0_xxyyy_yyzz[k] = -g_y_0_xyyy_yyzz[k] * cd_x[k] + g_y_0_xyyy_xyyzz[k];

                g_y_0_xxyyy_yzzz[k] = -g_y_0_xyyy_yzzz[k] * cd_x[k] + g_y_0_xyyy_xyzzz[k];

                g_y_0_xxyyy_zzzz[k] = -g_y_0_xyyy_zzzz[k] * cd_x[k] + g_y_0_xyyy_xzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 105);

            auto g_y_0_xxyyz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 106);

            auto g_y_0_xxyyz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 107);

            auto g_y_0_xxyyz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 108);

            auto g_y_0_xxyyz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 109);

            auto g_y_0_xxyyz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 110);

            auto g_y_0_xxyyz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 111);

            auto g_y_0_xxyyz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 112);

            auto g_y_0_xxyyz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 113);

            auto g_y_0_xxyyz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 114);

            auto g_y_0_xxyyz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 115);

            auto g_y_0_xxyyz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 116);

            auto g_y_0_xxyyz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 117);

            auto g_y_0_xxyyz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 118);

            auto g_y_0_xxyyz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyz_xxxx, g_y_0_xxyyz_xxxy, g_y_0_xxyyz_xxxz, g_y_0_xxyyz_xxyy, g_y_0_xxyyz_xxyz, g_y_0_xxyyz_xxzz, g_y_0_xxyyz_xyyy, g_y_0_xxyyz_xyyz, g_y_0_xxyyz_xyzz, g_y_0_xxyyz_xzzz, g_y_0_xxyyz_yyyy, g_y_0_xxyyz_yyyz, g_y_0_xxyyz_yyzz, g_y_0_xxyyz_yzzz, g_y_0_xxyyz_zzzz, g_y_0_xyyz_xxxx, g_y_0_xyyz_xxxxx, g_y_0_xyyz_xxxxy, g_y_0_xyyz_xxxxz, g_y_0_xyyz_xxxy, g_y_0_xyyz_xxxyy, g_y_0_xyyz_xxxyz, g_y_0_xyyz_xxxz, g_y_0_xyyz_xxxzz, g_y_0_xyyz_xxyy, g_y_0_xyyz_xxyyy, g_y_0_xyyz_xxyyz, g_y_0_xyyz_xxyz, g_y_0_xyyz_xxyzz, g_y_0_xyyz_xxzz, g_y_0_xyyz_xxzzz, g_y_0_xyyz_xyyy, g_y_0_xyyz_xyyyy, g_y_0_xyyz_xyyyz, g_y_0_xyyz_xyyz, g_y_0_xyyz_xyyzz, g_y_0_xyyz_xyzz, g_y_0_xyyz_xyzzz, g_y_0_xyyz_xzzz, g_y_0_xyyz_xzzzz, g_y_0_xyyz_yyyy, g_y_0_xyyz_yyyz, g_y_0_xyyz_yyzz, g_y_0_xyyz_yzzz, g_y_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_xxxx[k] = -g_y_0_xyyz_xxxx[k] * cd_x[k] + g_y_0_xyyz_xxxxx[k];

                g_y_0_xxyyz_xxxy[k] = -g_y_0_xyyz_xxxy[k] * cd_x[k] + g_y_0_xyyz_xxxxy[k];

                g_y_0_xxyyz_xxxz[k] = -g_y_0_xyyz_xxxz[k] * cd_x[k] + g_y_0_xyyz_xxxxz[k];

                g_y_0_xxyyz_xxyy[k] = -g_y_0_xyyz_xxyy[k] * cd_x[k] + g_y_0_xyyz_xxxyy[k];

                g_y_0_xxyyz_xxyz[k] = -g_y_0_xyyz_xxyz[k] * cd_x[k] + g_y_0_xyyz_xxxyz[k];

                g_y_0_xxyyz_xxzz[k] = -g_y_0_xyyz_xxzz[k] * cd_x[k] + g_y_0_xyyz_xxxzz[k];

                g_y_0_xxyyz_xyyy[k] = -g_y_0_xyyz_xyyy[k] * cd_x[k] + g_y_0_xyyz_xxyyy[k];

                g_y_0_xxyyz_xyyz[k] = -g_y_0_xyyz_xyyz[k] * cd_x[k] + g_y_0_xyyz_xxyyz[k];

                g_y_0_xxyyz_xyzz[k] = -g_y_0_xyyz_xyzz[k] * cd_x[k] + g_y_0_xyyz_xxyzz[k];

                g_y_0_xxyyz_xzzz[k] = -g_y_0_xyyz_xzzz[k] * cd_x[k] + g_y_0_xyyz_xxzzz[k];

                g_y_0_xxyyz_yyyy[k] = -g_y_0_xyyz_yyyy[k] * cd_x[k] + g_y_0_xyyz_xyyyy[k];

                g_y_0_xxyyz_yyyz[k] = -g_y_0_xyyz_yyyz[k] * cd_x[k] + g_y_0_xyyz_xyyyz[k];

                g_y_0_xxyyz_yyzz[k] = -g_y_0_xyyz_yyzz[k] * cd_x[k] + g_y_0_xyyz_xyyzz[k];

                g_y_0_xxyyz_yzzz[k] = -g_y_0_xyyz_yzzz[k] * cd_x[k] + g_y_0_xyyz_xyzzz[k];

                g_y_0_xxyyz_zzzz[k] = -g_y_0_xyyz_zzzz[k] * cd_x[k] + g_y_0_xyyz_xzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 120);

            auto g_y_0_xxyzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 121);

            auto g_y_0_xxyzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 122);

            auto g_y_0_xxyzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 123);

            auto g_y_0_xxyzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 124);

            auto g_y_0_xxyzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 125);

            auto g_y_0_xxyzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 126);

            auto g_y_0_xxyzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 127);

            auto g_y_0_xxyzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 128);

            auto g_y_0_xxyzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 129);

            auto g_y_0_xxyzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 130);

            auto g_y_0_xxyzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 131);

            auto g_y_0_xxyzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 132);

            auto g_y_0_xxyzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 133);

            auto g_y_0_xxyzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 134);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzz_xxxx, g_y_0_xxyzz_xxxy, g_y_0_xxyzz_xxxz, g_y_0_xxyzz_xxyy, g_y_0_xxyzz_xxyz, g_y_0_xxyzz_xxzz, g_y_0_xxyzz_xyyy, g_y_0_xxyzz_xyyz, g_y_0_xxyzz_xyzz, g_y_0_xxyzz_xzzz, g_y_0_xxyzz_yyyy, g_y_0_xxyzz_yyyz, g_y_0_xxyzz_yyzz, g_y_0_xxyzz_yzzz, g_y_0_xxyzz_zzzz, g_y_0_xyzz_xxxx, g_y_0_xyzz_xxxxx, g_y_0_xyzz_xxxxy, g_y_0_xyzz_xxxxz, g_y_0_xyzz_xxxy, g_y_0_xyzz_xxxyy, g_y_0_xyzz_xxxyz, g_y_0_xyzz_xxxz, g_y_0_xyzz_xxxzz, g_y_0_xyzz_xxyy, g_y_0_xyzz_xxyyy, g_y_0_xyzz_xxyyz, g_y_0_xyzz_xxyz, g_y_0_xyzz_xxyzz, g_y_0_xyzz_xxzz, g_y_0_xyzz_xxzzz, g_y_0_xyzz_xyyy, g_y_0_xyzz_xyyyy, g_y_0_xyzz_xyyyz, g_y_0_xyzz_xyyz, g_y_0_xyzz_xyyzz, g_y_0_xyzz_xyzz, g_y_0_xyzz_xyzzz, g_y_0_xyzz_xzzz, g_y_0_xyzz_xzzzz, g_y_0_xyzz_yyyy, g_y_0_xyzz_yyyz, g_y_0_xyzz_yyzz, g_y_0_xyzz_yzzz, g_y_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_xxxx[k] = -g_y_0_xyzz_xxxx[k] * cd_x[k] + g_y_0_xyzz_xxxxx[k];

                g_y_0_xxyzz_xxxy[k] = -g_y_0_xyzz_xxxy[k] * cd_x[k] + g_y_0_xyzz_xxxxy[k];

                g_y_0_xxyzz_xxxz[k] = -g_y_0_xyzz_xxxz[k] * cd_x[k] + g_y_0_xyzz_xxxxz[k];

                g_y_0_xxyzz_xxyy[k] = -g_y_0_xyzz_xxyy[k] * cd_x[k] + g_y_0_xyzz_xxxyy[k];

                g_y_0_xxyzz_xxyz[k] = -g_y_0_xyzz_xxyz[k] * cd_x[k] + g_y_0_xyzz_xxxyz[k];

                g_y_0_xxyzz_xxzz[k] = -g_y_0_xyzz_xxzz[k] * cd_x[k] + g_y_0_xyzz_xxxzz[k];

                g_y_0_xxyzz_xyyy[k] = -g_y_0_xyzz_xyyy[k] * cd_x[k] + g_y_0_xyzz_xxyyy[k];

                g_y_0_xxyzz_xyyz[k] = -g_y_0_xyzz_xyyz[k] * cd_x[k] + g_y_0_xyzz_xxyyz[k];

                g_y_0_xxyzz_xyzz[k] = -g_y_0_xyzz_xyzz[k] * cd_x[k] + g_y_0_xyzz_xxyzz[k];

                g_y_0_xxyzz_xzzz[k] = -g_y_0_xyzz_xzzz[k] * cd_x[k] + g_y_0_xyzz_xxzzz[k];

                g_y_0_xxyzz_yyyy[k] = -g_y_0_xyzz_yyyy[k] * cd_x[k] + g_y_0_xyzz_xyyyy[k];

                g_y_0_xxyzz_yyyz[k] = -g_y_0_xyzz_yyyz[k] * cd_x[k] + g_y_0_xyzz_xyyyz[k];

                g_y_0_xxyzz_yyzz[k] = -g_y_0_xyzz_yyzz[k] * cd_x[k] + g_y_0_xyzz_xyyzz[k];

                g_y_0_xxyzz_yzzz[k] = -g_y_0_xyzz_yzzz[k] * cd_x[k] + g_y_0_xyzz_xyzzz[k];

                g_y_0_xxyzz_zzzz[k] = -g_y_0_xyzz_zzzz[k] * cd_x[k] + g_y_0_xyzz_xzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 135);

            auto g_y_0_xxzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 136);

            auto g_y_0_xxzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 137);

            auto g_y_0_xxzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 138);

            auto g_y_0_xxzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 139);

            auto g_y_0_xxzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 140);

            auto g_y_0_xxzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 141);

            auto g_y_0_xxzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 142);

            auto g_y_0_xxzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 143);

            auto g_y_0_xxzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 144);

            auto g_y_0_xxzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 145);

            auto g_y_0_xxzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 146);

            auto g_y_0_xxzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 147);

            auto g_y_0_xxzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 148);

            auto g_y_0_xxzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzz_xxxx, g_y_0_xxzzz_xxxy, g_y_0_xxzzz_xxxz, g_y_0_xxzzz_xxyy, g_y_0_xxzzz_xxyz, g_y_0_xxzzz_xxzz, g_y_0_xxzzz_xyyy, g_y_0_xxzzz_xyyz, g_y_0_xxzzz_xyzz, g_y_0_xxzzz_xzzz, g_y_0_xxzzz_yyyy, g_y_0_xxzzz_yyyz, g_y_0_xxzzz_yyzz, g_y_0_xxzzz_yzzz, g_y_0_xxzzz_zzzz, g_y_0_xzzz_xxxx, g_y_0_xzzz_xxxxx, g_y_0_xzzz_xxxxy, g_y_0_xzzz_xxxxz, g_y_0_xzzz_xxxy, g_y_0_xzzz_xxxyy, g_y_0_xzzz_xxxyz, g_y_0_xzzz_xxxz, g_y_0_xzzz_xxxzz, g_y_0_xzzz_xxyy, g_y_0_xzzz_xxyyy, g_y_0_xzzz_xxyyz, g_y_0_xzzz_xxyz, g_y_0_xzzz_xxyzz, g_y_0_xzzz_xxzz, g_y_0_xzzz_xxzzz, g_y_0_xzzz_xyyy, g_y_0_xzzz_xyyyy, g_y_0_xzzz_xyyyz, g_y_0_xzzz_xyyz, g_y_0_xzzz_xyyzz, g_y_0_xzzz_xyzz, g_y_0_xzzz_xyzzz, g_y_0_xzzz_xzzz, g_y_0_xzzz_xzzzz, g_y_0_xzzz_yyyy, g_y_0_xzzz_yyyz, g_y_0_xzzz_yyzz, g_y_0_xzzz_yzzz, g_y_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_xxxx[k] = -g_y_0_xzzz_xxxx[k] * cd_x[k] + g_y_0_xzzz_xxxxx[k];

                g_y_0_xxzzz_xxxy[k] = -g_y_0_xzzz_xxxy[k] * cd_x[k] + g_y_0_xzzz_xxxxy[k];

                g_y_0_xxzzz_xxxz[k] = -g_y_0_xzzz_xxxz[k] * cd_x[k] + g_y_0_xzzz_xxxxz[k];

                g_y_0_xxzzz_xxyy[k] = -g_y_0_xzzz_xxyy[k] * cd_x[k] + g_y_0_xzzz_xxxyy[k];

                g_y_0_xxzzz_xxyz[k] = -g_y_0_xzzz_xxyz[k] * cd_x[k] + g_y_0_xzzz_xxxyz[k];

                g_y_0_xxzzz_xxzz[k] = -g_y_0_xzzz_xxzz[k] * cd_x[k] + g_y_0_xzzz_xxxzz[k];

                g_y_0_xxzzz_xyyy[k] = -g_y_0_xzzz_xyyy[k] * cd_x[k] + g_y_0_xzzz_xxyyy[k];

                g_y_0_xxzzz_xyyz[k] = -g_y_0_xzzz_xyyz[k] * cd_x[k] + g_y_0_xzzz_xxyyz[k];

                g_y_0_xxzzz_xyzz[k] = -g_y_0_xzzz_xyzz[k] * cd_x[k] + g_y_0_xzzz_xxyzz[k];

                g_y_0_xxzzz_xzzz[k] = -g_y_0_xzzz_xzzz[k] * cd_x[k] + g_y_0_xzzz_xxzzz[k];

                g_y_0_xxzzz_yyyy[k] = -g_y_0_xzzz_yyyy[k] * cd_x[k] + g_y_0_xzzz_xyyyy[k];

                g_y_0_xxzzz_yyyz[k] = -g_y_0_xzzz_yyyz[k] * cd_x[k] + g_y_0_xzzz_xyyyz[k];

                g_y_0_xxzzz_yyzz[k] = -g_y_0_xzzz_yyzz[k] * cd_x[k] + g_y_0_xzzz_xyyzz[k];

                g_y_0_xxzzz_yzzz[k] = -g_y_0_xzzz_yzzz[k] * cd_x[k] + g_y_0_xzzz_xyzzz[k];

                g_y_0_xxzzz_zzzz[k] = -g_y_0_xzzz_zzzz[k] * cd_x[k] + g_y_0_xzzz_xzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyy_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 150);

            auto g_y_0_xyyyy_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 151);

            auto g_y_0_xyyyy_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 152);

            auto g_y_0_xyyyy_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 153);

            auto g_y_0_xyyyy_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 154);

            auto g_y_0_xyyyy_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 155);

            auto g_y_0_xyyyy_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 156);

            auto g_y_0_xyyyy_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 157);

            auto g_y_0_xyyyy_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 158);

            auto g_y_0_xyyyy_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 159);

            auto g_y_0_xyyyy_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 160);

            auto g_y_0_xyyyy_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 161);

            auto g_y_0_xyyyy_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 162);

            auto g_y_0_xyyyy_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 163);

            auto g_y_0_xyyyy_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 164);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyy_xxxx, g_y_0_xyyyy_xxxy, g_y_0_xyyyy_xxxz, g_y_0_xyyyy_xxyy, g_y_0_xyyyy_xxyz, g_y_0_xyyyy_xxzz, g_y_0_xyyyy_xyyy, g_y_0_xyyyy_xyyz, g_y_0_xyyyy_xyzz, g_y_0_xyyyy_xzzz, g_y_0_xyyyy_yyyy, g_y_0_xyyyy_yyyz, g_y_0_xyyyy_yyzz, g_y_0_xyyyy_yzzz, g_y_0_xyyyy_zzzz, g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxxx, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_xxxx[k] = -g_y_0_yyyy_xxxx[k] * cd_x[k] + g_y_0_yyyy_xxxxx[k];

                g_y_0_xyyyy_xxxy[k] = -g_y_0_yyyy_xxxy[k] * cd_x[k] + g_y_0_yyyy_xxxxy[k];

                g_y_0_xyyyy_xxxz[k] = -g_y_0_yyyy_xxxz[k] * cd_x[k] + g_y_0_yyyy_xxxxz[k];

                g_y_0_xyyyy_xxyy[k] = -g_y_0_yyyy_xxyy[k] * cd_x[k] + g_y_0_yyyy_xxxyy[k];

                g_y_0_xyyyy_xxyz[k] = -g_y_0_yyyy_xxyz[k] * cd_x[k] + g_y_0_yyyy_xxxyz[k];

                g_y_0_xyyyy_xxzz[k] = -g_y_0_yyyy_xxzz[k] * cd_x[k] + g_y_0_yyyy_xxxzz[k];

                g_y_0_xyyyy_xyyy[k] = -g_y_0_yyyy_xyyy[k] * cd_x[k] + g_y_0_yyyy_xxyyy[k];

                g_y_0_xyyyy_xyyz[k] = -g_y_0_yyyy_xyyz[k] * cd_x[k] + g_y_0_yyyy_xxyyz[k];

                g_y_0_xyyyy_xyzz[k] = -g_y_0_yyyy_xyzz[k] * cd_x[k] + g_y_0_yyyy_xxyzz[k];

                g_y_0_xyyyy_xzzz[k] = -g_y_0_yyyy_xzzz[k] * cd_x[k] + g_y_0_yyyy_xxzzz[k];

                g_y_0_xyyyy_yyyy[k] = -g_y_0_yyyy_yyyy[k] * cd_x[k] + g_y_0_yyyy_xyyyy[k];

                g_y_0_xyyyy_yyyz[k] = -g_y_0_yyyy_yyyz[k] * cd_x[k] + g_y_0_yyyy_xyyyz[k];

                g_y_0_xyyyy_yyzz[k] = -g_y_0_yyyy_yyzz[k] * cd_x[k] + g_y_0_yyyy_xyyzz[k];

                g_y_0_xyyyy_yzzz[k] = -g_y_0_yyyy_yzzz[k] * cd_x[k] + g_y_0_yyyy_xyzzz[k];

                g_y_0_xyyyy_zzzz[k] = -g_y_0_yyyy_zzzz[k] * cd_x[k] + g_y_0_yyyy_xzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 165);

            auto g_y_0_xyyyz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 166);

            auto g_y_0_xyyyz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 167);

            auto g_y_0_xyyyz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 168);

            auto g_y_0_xyyyz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 169);

            auto g_y_0_xyyyz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 170);

            auto g_y_0_xyyyz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 171);

            auto g_y_0_xyyyz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 172);

            auto g_y_0_xyyyz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 173);

            auto g_y_0_xyyyz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 174);

            auto g_y_0_xyyyz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 175);

            auto g_y_0_xyyyz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 176);

            auto g_y_0_xyyyz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 177);

            auto g_y_0_xyyyz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 178);

            auto g_y_0_xyyyz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyz_xxxx, g_y_0_xyyyz_xxxy, g_y_0_xyyyz_xxxz, g_y_0_xyyyz_xxyy, g_y_0_xyyyz_xxyz, g_y_0_xyyyz_xxzz, g_y_0_xyyyz_xyyy, g_y_0_xyyyz_xyyz, g_y_0_xyyyz_xyzz, g_y_0_xyyyz_xzzz, g_y_0_xyyyz_yyyy, g_y_0_xyyyz_yyyz, g_y_0_xyyyz_yyzz, g_y_0_xyyyz_yzzz, g_y_0_xyyyz_zzzz, g_y_0_yyyz_xxxx, g_y_0_yyyz_xxxxx, g_y_0_yyyz_xxxxy, g_y_0_yyyz_xxxxz, g_y_0_yyyz_xxxy, g_y_0_yyyz_xxxyy, g_y_0_yyyz_xxxyz, g_y_0_yyyz_xxxz, g_y_0_yyyz_xxxzz, g_y_0_yyyz_xxyy, g_y_0_yyyz_xxyyy, g_y_0_yyyz_xxyyz, g_y_0_yyyz_xxyz, g_y_0_yyyz_xxyzz, g_y_0_yyyz_xxzz, g_y_0_yyyz_xxzzz, g_y_0_yyyz_xyyy, g_y_0_yyyz_xyyyy, g_y_0_yyyz_xyyyz, g_y_0_yyyz_xyyz, g_y_0_yyyz_xyyzz, g_y_0_yyyz_xyzz, g_y_0_yyyz_xyzzz, g_y_0_yyyz_xzzz, g_y_0_yyyz_xzzzz, g_y_0_yyyz_yyyy, g_y_0_yyyz_yyyz, g_y_0_yyyz_yyzz, g_y_0_yyyz_yzzz, g_y_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_xxxx[k] = -g_y_0_yyyz_xxxx[k] * cd_x[k] + g_y_0_yyyz_xxxxx[k];

                g_y_0_xyyyz_xxxy[k] = -g_y_0_yyyz_xxxy[k] * cd_x[k] + g_y_0_yyyz_xxxxy[k];

                g_y_0_xyyyz_xxxz[k] = -g_y_0_yyyz_xxxz[k] * cd_x[k] + g_y_0_yyyz_xxxxz[k];

                g_y_0_xyyyz_xxyy[k] = -g_y_0_yyyz_xxyy[k] * cd_x[k] + g_y_0_yyyz_xxxyy[k];

                g_y_0_xyyyz_xxyz[k] = -g_y_0_yyyz_xxyz[k] * cd_x[k] + g_y_0_yyyz_xxxyz[k];

                g_y_0_xyyyz_xxzz[k] = -g_y_0_yyyz_xxzz[k] * cd_x[k] + g_y_0_yyyz_xxxzz[k];

                g_y_0_xyyyz_xyyy[k] = -g_y_0_yyyz_xyyy[k] * cd_x[k] + g_y_0_yyyz_xxyyy[k];

                g_y_0_xyyyz_xyyz[k] = -g_y_0_yyyz_xyyz[k] * cd_x[k] + g_y_0_yyyz_xxyyz[k];

                g_y_0_xyyyz_xyzz[k] = -g_y_0_yyyz_xyzz[k] * cd_x[k] + g_y_0_yyyz_xxyzz[k];

                g_y_0_xyyyz_xzzz[k] = -g_y_0_yyyz_xzzz[k] * cd_x[k] + g_y_0_yyyz_xxzzz[k];

                g_y_0_xyyyz_yyyy[k] = -g_y_0_yyyz_yyyy[k] * cd_x[k] + g_y_0_yyyz_xyyyy[k];

                g_y_0_xyyyz_yyyz[k] = -g_y_0_yyyz_yyyz[k] * cd_x[k] + g_y_0_yyyz_xyyyz[k];

                g_y_0_xyyyz_yyzz[k] = -g_y_0_yyyz_yyzz[k] * cd_x[k] + g_y_0_yyyz_xyyzz[k];

                g_y_0_xyyyz_yzzz[k] = -g_y_0_yyyz_yzzz[k] * cd_x[k] + g_y_0_yyyz_xyzzz[k];

                g_y_0_xyyyz_zzzz[k] = -g_y_0_yyyz_zzzz[k] * cd_x[k] + g_y_0_yyyz_xzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 180);

            auto g_y_0_xyyzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 181);

            auto g_y_0_xyyzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 182);

            auto g_y_0_xyyzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 183);

            auto g_y_0_xyyzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 184);

            auto g_y_0_xyyzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 185);

            auto g_y_0_xyyzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 186);

            auto g_y_0_xyyzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 187);

            auto g_y_0_xyyzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 188);

            auto g_y_0_xyyzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 189);

            auto g_y_0_xyyzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 190);

            auto g_y_0_xyyzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 191);

            auto g_y_0_xyyzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 192);

            auto g_y_0_xyyzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 193);

            auto g_y_0_xyyzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 194);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzz_xxxx, g_y_0_xyyzz_xxxy, g_y_0_xyyzz_xxxz, g_y_0_xyyzz_xxyy, g_y_0_xyyzz_xxyz, g_y_0_xyyzz_xxzz, g_y_0_xyyzz_xyyy, g_y_0_xyyzz_xyyz, g_y_0_xyyzz_xyzz, g_y_0_xyyzz_xzzz, g_y_0_xyyzz_yyyy, g_y_0_xyyzz_yyyz, g_y_0_xyyzz_yyzz, g_y_0_xyyzz_yzzz, g_y_0_xyyzz_zzzz, g_y_0_yyzz_xxxx, g_y_0_yyzz_xxxxx, g_y_0_yyzz_xxxxy, g_y_0_yyzz_xxxxz, g_y_0_yyzz_xxxy, g_y_0_yyzz_xxxyy, g_y_0_yyzz_xxxyz, g_y_0_yyzz_xxxz, g_y_0_yyzz_xxxzz, g_y_0_yyzz_xxyy, g_y_0_yyzz_xxyyy, g_y_0_yyzz_xxyyz, g_y_0_yyzz_xxyz, g_y_0_yyzz_xxyzz, g_y_0_yyzz_xxzz, g_y_0_yyzz_xxzzz, g_y_0_yyzz_xyyy, g_y_0_yyzz_xyyyy, g_y_0_yyzz_xyyyz, g_y_0_yyzz_xyyz, g_y_0_yyzz_xyyzz, g_y_0_yyzz_xyzz, g_y_0_yyzz_xyzzz, g_y_0_yyzz_xzzz, g_y_0_yyzz_xzzzz, g_y_0_yyzz_yyyy, g_y_0_yyzz_yyyz, g_y_0_yyzz_yyzz, g_y_0_yyzz_yzzz, g_y_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_xxxx[k] = -g_y_0_yyzz_xxxx[k] * cd_x[k] + g_y_0_yyzz_xxxxx[k];

                g_y_0_xyyzz_xxxy[k] = -g_y_0_yyzz_xxxy[k] * cd_x[k] + g_y_0_yyzz_xxxxy[k];

                g_y_0_xyyzz_xxxz[k] = -g_y_0_yyzz_xxxz[k] * cd_x[k] + g_y_0_yyzz_xxxxz[k];

                g_y_0_xyyzz_xxyy[k] = -g_y_0_yyzz_xxyy[k] * cd_x[k] + g_y_0_yyzz_xxxyy[k];

                g_y_0_xyyzz_xxyz[k] = -g_y_0_yyzz_xxyz[k] * cd_x[k] + g_y_0_yyzz_xxxyz[k];

                g_y_0_xyyzz_xxzz[k] = -g_y_0_yyzz_xxzz[k] * cd_x[k] + g_y_0_yyzz_xxxzz[k];

                g_y_0_xyyzz_xyyy[k] = -g_y_0_yyzz_xyyy[k] * cd_x[k] + g_y_0_yyzz_xxyyy[k];

                g_y_0_xyyzz_xyyz[k] = -g_y_0_yyzz_xyyz[k] * cd_x[k] + g_y_0_yyzz_xxyyz[k];

                g_y_0_xyyzz_xyzz[k] = -g_y_0_yyzz_xyzz[k] * cd_x[k] + g_y_0_yyzz_xxyzz[k];

                g_y_0_xyyzz_xzzz[k] = -g_y_0_yyzz_xzzz[k] * cd_x[k] + g_y_0_yyzz_xxzzz[k];

                g_y_0_xyyzz_yyyy[k] = -g_y_0_yyzz_yyyy[k] * cd_x[k] + g_y_0_yyzz_xyyyy[k];

                g_y_0_xyyzz_yyyz[k] = -g_y_0_yyzz_yyyz[k] * cd_x[k] + g_y_0_yyzz_xyyyz[k];

                g_y_0_xyyzz_yyzz[k] = -g_y_0_yyzz_yyzz[k] * cd_x[k] + g_y_0_yyzz_xyyzz[k];

                g_y_0_xyyzz_yzzz[k] = -g_y_0_yyzz_yzzz[k] * cd_x[k] + g_y_0_yyzz_xyzzz[k];

                g_y_0_xyyzz_zzzz[k] = -g_y_0_yyzz_zzzz[k] * cd_x[k] + g_y_0_yyzz_xzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 195);

            auto g_y_0_xyzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 196);

            auto g_y_0_xyzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 197);

            auto g_y_0_xyzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 198);

            auto g_y_0_xyzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 199);

            auto g_y_0_xyzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 200);

            auto g_y_0_xyzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 201);

            auto g_y_0_xyzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 202);

            auto g_y_0_xyzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 203);

            auto g_y_0_xyzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 204);

            auto g_y_0_xyzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 205);

            auto g_y_0_xyzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 206);

            auto g_y_0_xyzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 207);

            auto g_y_0_xyzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 208);

            auto g_y_0_xyzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzz_xxxx, g_y_0_xyzzz_xxxy, g_y_0_xyzzz_xxxz, g_y_0_xyzzz_xxyy, g_y_0_xyzzz_xxyz, g_y_0_xyzzz_xxzz, g_y_0_xyzzz_xyyy, g_y_0_xyzzz_xyyz, g_y_0_xyzzz_xyzz, g_y_0_xyzzz_xzzz, g_y_0_xyzzz_yyyy, g_y_0_xyzzz_yyyz, g_y_0_xyzzz_yyzz, g_y_0_xyzzz_yzzz, g_y_0_xyzzz_zzzz, g_y_0_yzzz_xxxx, g_y_0_yzzz_xxxxx, g_y_0_yzzz_xxxxy, g_y_0_yzzz_xxxxz, g_y_0_yzzz_xxxy, g_y_0_yzzz_xxxyy, g_y_0_yzzz_xxxyz, g_y_0_yzzz_xxxz, g_y_0_yzzz_xxxzz, g_y_0_yzzz_xxyy, g_y_0_yzzz_xxyyy, g_y_0_yzzz_xxyyz, g_y_0_yzzz_xxyz, g_y_0_yzzz_xxyzz, g_y_0_yzzz_xxzz, g_y_0_yzzz_xxzzz, g_y_0_yzzz_xyyy, g_y_0_yzzz_xyyyy, g_y_0_yzzz_xyyyz, g_y_0_yzzz_xyyz, g_y_0_yzzz_xyyzz, g_y_0_yzzz_xyzz, g_y_0_yzzz_xyzzz, g_y_0_yzzz_xzzz, g_y_0_yzzz_xzzzz, g_y_0_yzzz_yyyy, g_y_0_yzzz_yyyz, g_y_0_yzzz_yyzz, g_y_0_yzzz_yzzz, g_y_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_xxxx[k] = -g_y_0_yzzz_xxxx[k] * cd_x[k] + g_y_0_yzzz_xxxxx[k];

                g_y_0_xyzzz_xxxy[k] = -g_y_0_yzzz_xxxy[k] * cd_x[k] + g_y_0_yzzz_xxxxy[k];

                g_y_0_xyzzz_xxxz[k] = -g_y_0_yzzz_xxxz[k] * cd_x[k] + g_y_0_yzzz_xxxxz[k];

                g_y_0_xyzzz_xxyy[k] = -g_y_0_yzzz_xxyy[k] * cd_x[k] + g_y_0_yzzz_xxxyy[k];

                g_y_0_xyzzz_xxyz[k] = -g_y_0_yzzz_xxyz[k] * cd_x[k] + g_y_0_yzzz_xxxyz[k];

                g_y_0_xyzzz_xxzz[k] = -g_y_0_yzzz_xxzz[k] * cd_x[k] + g_y_0_yzzz_xxxzz[k];

                g_y_0_xyzzz_xyyy[k] = -g_y_0_yzzz_xyyy[k] * cd_x[k] + g_y_0_yzzz_xxyyy[k];

                g_y_0_xyzzz_xyyz[k] = -g_y_0_yzzz_xyyz[k] * cd_x[k] + g_y_0_yzzz_xxyyz[k];

                g_y_0_xyzzz_xyzz[k] = -g_y_0_yzzz_xyzz[k] * cd_x[k] + g_y_0_yzzz_xxyzz[k];

                g_y_0_xyzzz_xzzz[k] = -g_y_0_yzzz_xzzz[k] * cd_x[k] + g_y_0_yzzz_xxzzz[k];

                g_y_0_xyzzz_yyyy[k] = -g_y_0_yzzz_yyyy[k] * cd_x[k] + g_y_0_yzzz_xyyyy[k];

                g_y_0_xyzzz_yyyz[k] = -g_y_0_yzzz_yyyz[k] * cd_x[k] + g_y_0_yzzz_xyyyz[k];

                g_y_0_xyzzz_yyzz[k] = -g_y_0_yzzz_yyzz[k] * cd_x[k] + g_y_0_yzzz_xyyzz[k];

                g_y_0_xyzzz_yzzz[k] = -g_y_0_yzzz_yzzz[k] * cd_x[k] + g_y_0_yzzz_xyzzz[k];

                g_y_0_xyzzz_zzzz[k] = -g_y_0_yzzz_zzzz[k] * cd_x[k] + g_y_0_yzzz_xzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 210);

            auto g_y_0_xzzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 211);

            auto g_y_0_xzzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 212);

            auto g_y_0_xzzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 213);

            auto g_y_0_xzzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 214);

            auto g_y_0_xzzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 215);

            auto g_y_0_xzzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 216);

            auto g_y_0_xzzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 217);

            auto g_y_0_xzzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 218);

            auto g_y_0_xzzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 219);

            auto g_y_0_xzzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 220);

            auto g_y_0_xzzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 221);

            auto g_y_0_xzzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 222);

            auto g_y_0_xzzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 223);

            auto g_y_0_xzzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 224);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzz_xxxx, g_y_0_xzzzz_xxxy, g_y_0_xzzzz_xxxz, g_y_0_xzzzz_xxyy, g_y_0_xzzzz_xxyz, g_y_0_xzzzz_xxzz, g_y_0_xzzzz_xyyy, g_y_0_xzzzz_xyyz, g_y_0_xzzzz_xyzz, g_y_0_xzzzz_xzzz, g_y_0_xzzzz_yyyy, g_y_0_xzzzz_yyyz, g_y_0_xzzzz_yyzz, g_y_0_xzzzz_yzzz, g_y_0_xzzzz_zzzz, g_y_0_zzzz_xxxx, g_y_0_zzzz_xxxxx, g_y_0_zzzz_xxxxy, g_y_0_zzzz_xxxxz, g_y_0_zzzz_xxxy, g_y_0_zzzz_xxxyy, g_y_0_zzzz_xxxyz, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxxzz, g_y_0_zzzz_xxyy, g_y_0_zzzz_xxyyy, g_y_0_zzzz_xxyyz, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxyzz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xxzzz, g_y_0_zzzz_xyyy, g_y_0_zzzz_xyyyy, g_y_0_zzzz_xyyyz, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyyzz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xyzzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_xzzzz, g_y_0_zzzz_yyyy, g_y_0_zzzz_yyyz, g_y_0_zzzz_yyzz, g_y_0_zzzz_yzzz, g_y_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_xxxx[k] = -g_y_0_zzzz_xxxx[k] * cd_x[k] + g_y_0_zzzz_xxxxx[k];

                g_y_0_xzzzz_xxxy[k] = -g_y_0_zzzz_xxxy[k] * cd_x[k] + g_y_0_zzzz_xxxxy[k];

                g_y_0_xzzzz_xxxz[k] = -g_y_0_zzzz_xxxz[k] * cd_x[k] + g_y_0_zzzz_xxxxz[k];

                g_y_0_xzzzz_xxyy[k] = -g_y_0_zzzz_xxyy[k] * cd_x[k] + g_y_0_zzzz_xxxyy[k];

                g_y_0_xzzzz_xxyz[k] = -g_y_0_zzzz_xxyz[k] * cd_x[k] + g_y_0_zzzz_xxxyz[k];

                g_y_0_xzzzz_xxzz[k] = -g_y_0_zzzz_xxzz[k] * cd_x[k] + g_y_0_zzzz_xxxzz[k];

                g_y_0_xzzzz_xyyy[k] = -g_y_0_zzzz_xyyy[k] * cd_x[k] + g_y_0_zzzz_xxyyy[k];

                g_y_0_xzzzz_xyyz[k] = -g_y_0_zzzz_xyyz[k] * cd_x[k] + g_y_0_zzzz_xxyyz[k];

                g_y_0_xzzzz_xyzz[k] = -g_y_0_zzzz_xyzz[k] * cd_x[k] + g_y_0_zzzz_xxyzz[k];

                g_y_0_xzzzz_xzzz[k] = -g_y_0_zzzz_xzzz[k] * cd_x[k] + g_y_0_zzzz_xxzzz[k];

                g_y_0_xzzzz_yyyy[k] = -g_y_0_zzzz_yyyy[k] * cd_x[k] + g_y_0_zzzz_xyyyy[k];

                g_y_0_xzzzz_yyyz[k] = -g_y_0_zzzz_yyyz[k] * cd_x[k] + g_y_0_zzzz_xyyyz[k];

                g_y_0_xzzzz_yyzz[k] = -g_y_0_zzzz_yyzz[k] * cd_x[k] + g_y_0_zzzz_xyyzz[k];

                g_y_0_xzzzz_yzzz[k] = -g_y_0_zzzz_yzzz[k] * cd_x[k] + g_y_0_zzzz_xyzzz[k];

                g_y_0_xzzzz_zzzz[k] = -g_y_0_zzzz_zzzz[k] * cd_x[k] + g_y_0_zzzz_xzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyy_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 225);

            auto g_y_0_yyyyy_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 226);

            auto g_y_0_yyyyy_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 227);

            auto g_y_0_yyyyy_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 228);

            auto g_y_0_yyyyy_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 229);

            auto g_y_0_yyyyy_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 230);

            auto g_y_0_yyyyy_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 231);

            auto g_y_0_yyyyy_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 232);

            auto g_y_0_yyyyy_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 233);

            auto g_y_0_yyyyy_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 234);

            auto g_y_0_yyyyy_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 235);

            auto g_y_0_yyyyy_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 236);

            auto g_y_0_yyyyy_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 237);

            auto g_y_0_yyyyy_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 238);

            auto g_y_0_yyyyy_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 239);

            #pragma omp simd aligned(cd_y, g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyyy, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_zzzz, g_y_0_yyyyy_xxxx, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_yyyy, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_zzzz, g_yyyy_xxxx, g_yyyy_xxxy, g_yyyy_xxxz, g_yyyy_xxyy, g_yyyy_xxyz, g_yyyy_xxzz, g_yyyy_xyyy, g_yyyy_xyyz, g_yyyy_xyzz, g_yyyy_xzzz, g_yyyy_yyyy, g_yyyy_yyyz, g_yyyy_yyzz, g_yyyy_yzzz, g_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_xxxx[k] = -g_yyyy_xxxx[k] - g_y_0_yyyy_xxxx[k] * cd_y[k] + g_y_0_yyyy_xxxxy[k];

                g_y_0_yyyyy_xxxy[k] = -g_yyyy_xxxy[k] - g_y_0_yyyy_xxxy[k] * cd_y[k] + g_y_0_yyyy_xxxyy[k];

                g_y_0_yyyyy_xxxz[k] = -g_yyyy_xxxz[k] - g_y_0_yyyy_xxxz[k] * cd_y[k] + g_y_0_yyyy_xxxyz[k];

                g_y_0_yyyyy_xxyy[k] = -g_yyyy_xxyy[k] - g_y_0_yyyy_xxyy[k] * cd_y[k] + g_y_0_yyyy_xxyyy[k];

                g_y_0_yyyyy_xxyz[k] = -g_yyyy_xxyz[k] - g_y_0_yyyy_xxyz[k] * cd_y[k] + g_y_0_yyyy_xxyyz[k];

                g_y_0_yyyyy_xxzz[k] = -g_yyyy_xxzz[k] - g_y_0_yyyy_xxzz[k] * cd_y[k] + g_y_0_yyyy_xxyzz[k];

                g_y_0_yyyyy_xyyy[k] = -g_yyyy_xyyy[k] - g_y_0_yyyy_xyyy[k] * cd_y[k] + g_y_0_yyyy_xyyyy[k];

                g_y_0_yyyyy_xyyz[k] = -g_yyyy_xyyz[k] - g_y_0_yyyy_xyyz[k] * cd_y[k] + g_y_0_yyyy_xyyyz[k];

                g_y_0_yyyyy_xyzz[k] = -g_yyyy_xyzz[k] - g_y_0_yyyy_xyzz[k] * cd_y[k] + g_y_0_yyyy_xyyzz[k];

                g_y_0_yyyyy_xzzz[k] = -g_yyyy_xzzz[k] - g_y_0_yyyy_xzzz[k] * cd_y[k] + g_y_0_yyyy_xyzzz[k];

                g_y_0_yyyyy_yyyy[k] = -g_yyyy_yyyy[k] - g_y_0_yyyy_yyyy[k] * cd_y[k] + g_y_0_yyyy_yyyyy[k];

                g_y_0_yyyyy_yyyz[k] = -g_yyyy_yyyz[k] - g_y_0_yyyy_yyyz[k] * cd_y[k] + g_y_0_yyyy_yyyyz[k];

                g_y_0_yyyyy_yyzz[k] = -g_yyyy_yyzz[k] - g_y_0_yyyy_yyzz[k] * cd_y[k] + g_y_0_yyyy_yyyzz[k];

                g_y_0_yyyyy_yzzz[k] = -g_yyyy_yzzz[k] - g_y_0_yyyy_yzzz[k] * cd_y[k] + g_y_0_yyyy_yyzzz[k];

                g_y_0_yyyyy_zzzz[k] = -g_yyyy_zzzz[k] - g_y_0_yyyy_zzzz[k] * cd_y[k] + g_y_0_yyyy_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 240);

            auto g_y_0_yyyyz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 241);

            auto g_y_0_yyyyz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 242);

            auto g_y_0_yyyyz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 243);

            auto g_y_0_yyyyz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 244);

            auto g_y_0_yyyyz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 245);

            auto g_y_0_yyyyz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 246);

            auto g_y_0_yyyyz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 247);

            auto g_y_0_yyyyz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 248);

            auto g_y_0_yyyyz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 249);

            auto g_y_0_yyyyz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 250);

            auto g_y_0_yyyyz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 251);

            auto g_y_0_yyyyz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 252);

            auto g_y_0_yyyyz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 253);

            auto g_y_0_yyyyz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 254);

            #pragma omp simd aligned(cd_z, g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_zzzz, g_y_0_yyyy_zzzzz, g_y_0_yyyyz_xxxx, g_y_0_yyyyz_xxxy, g_y_0_yyyyz_xxxz, g_y_0_yyyyz_xxyy, g_y_0_yyyyz_xxyz, g_y_0_yyyyz_xxzz, g_y_0_yyyyz_xyyy, g_y_0_yyyyz_xyyz, g_y_0_yyyyz_xyzz, g_y_0_yyyyz_xzzz, g_y_0_yyyyz_yyyy, g_y_0_yyyyz_yyyz, g_y_0_yyyyz_yyzz, g_y_0_yyyyz_yzzz, g_y_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_xxxx[k] = -g_y_0_yyyy_xxxx[k] * cd_z[k] + g_y_0_yyyy_xxxxz[k];

                g_y_0_yyyyz_xxxy[k] = -g_y_0_yyyy_xxxy[k] * cd_z[k] + g_y_0_yyyy_xxxyz[k];

                g_y_0_yyyyz_xxxz[k] = -g_y_0_yyyy_xxxz[k] * cd_z[k] + g_y_0_yyyy_xxxzz[k];

                g_y_0_yyyyz_xxyy[k] = -g_y_0_yyyy_xxyy[k] * cd_z[k] + g_y_0_yyyy_xxyyz[k];

                g_y_0_yyyyz_xxyz[k] = -g_y_0_yyyy_xxyz[k] * cd_z[k] + g_y_0_yyyy_xxyzz[k];

                g_y_0_yyyyz_xxzz[k] = -g_y_0_yyyy_xxzz[k] * cd_z[k] + g_y_0_yyyy_xxzzz[k];

                g_y_0_yyyyz_xyyy[k] = -g_y_0_yyyy_xyyy[k] * cd_z[k] + g_y_0_yyyy_xyyyz[k];

                g_y_0_yyyyz_xyyz[k] = -g_y_0_yyyy_xyyz[k] * cd_z[k] + g_y_0_yyyy_xyyzz[k];

                g_y_0_yyyyz_xyzz[k] = -g_y_0_yyyy_xyzz[k] * cd_z[k] + g_y_0_yyyy_xyzzz[k];

                g_y_0_yyyyz_xzzz[k] = -g_y_0_yyyy_xzzz[k] * cd_z[k] + g_y_0_yyyy_xzzzz[k];

                g_y_0_yyyyz_yyyy[k] = -g_y_0_yyyy_yyyy[k] * cd_z[k] + g_y_0_yyyy_yyyyz[k];

                g_y_0_yyyyz_yyyz[k] = -g_y_0_yyyy_yyyz[k] * cd_z[k] + g_y_0_yyyy_yyyzz[k];

                g_y_0_yyyyz_yyzz[k] = -g_y_0_yyyy_yyzz[k] * cd_z[k] + g_y_0_yyyy_yyzzz[k];

                g_y_0_yyyyz_yzzz[k] = -g_y_0_yyyy_yzzz[k] * cd_z[k] + g_y_0_yyyy_yzzzz[k];

                g_y_0_yyyyz_zzzz[k] = -g_y_0_yyyy_zzzz[k] * cd_z[k] + g_y_0_yyyy_zzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 255);

            auto g_y_0_yyyzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 256);

            auto g_y_0_yyyzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 257);

            auto g_y_0_yyyzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 258);

            auto g_y_0_yyyzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 259);

            auto g_y_0_yyyzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 260);

            auto g_y_0_yyyzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 261);

            auto g_y_0_yyyzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 262);

            auto g_y_0_yyyzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 263);

            auto g_y_0_yyyzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 264);

            auto g_y_0_yyyzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 265);

            auto g_y_0_yyyzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 266);

            auto g_y_0_yyyzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 267);

            auto g_y_0_yyyzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 268);

            auto g_y_0_yyyzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 269);

            #pragma omp simd aligned(cd_z, g_y_0_yyyz_xxxx, g_y_0_yyyz_xxxxz, g_y_0_yyyz_xxxy, g_y_0_yyyz_xxxyz, g_y_0_yyyz_xxxz, g_y_0_yyyz_xxxzz, g_y_0_yyyz_xxyy, g_y_0_yyyz_xxyyz, g_y_0_yyyz_xxyz, g_y_0_yyyz_xxyzz, g_y_0_yyyz_xxzz, g_y_0_yyyz_xxzzz, g_y_0_yyyz_xyyy, g_y_0_yyyz_xyyyz, g_y_0_yyyz_xyyz, g_y_0_yyyz_xyyzz, g_y_0_yyyz_xyzz, g_y_0_yyyz_xyzzz, g_y_0_yyyz_xzzz, g_y_0_yyyz_xzzzz, g_y_0_yyyz_yyyy, g_y_0_yyyz_yyyyz, g_y_0_yyyz_yyyz, g_y_0_yyyz_yyyzz, g_y_0_yyyz_yyzz, g_y_0_yyyz_yyzzz, g_y_0_yyyz_yzzz, g_y_0_yyyz_yzzzz, g_y_0_yyyz_zzzz, g_y_0_yyyz_zzzzz, g_y_0_yyyzz_xxxx, g_y_0_yyyzz_xxxy, g_y_0_yyyzz_xxxz, g_y_0_yyyzz_xxyy, g_y_0_yyyzz_xxyz, g_y_0_yyyzz_xxzz, g_y_0_yyyzz_xyyy, g_y_0_yyyzz_xyyz, g_y_0_yyyzz_xyzz, g_y_0_yyyzz_xzzz, g_y_0_yyyzz_yyyy, g_y_0_yyyzz_yyyz, g_y_0_yyyzz_yyzz, g_y_0_yyyzz_yzzz, g_y_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_xxxx[k] = -g_y_0_yyyz_xxxx[k] * cd_z[k] + g_y_0_yyyz_xxxxz[k];

                g_y_0_yyyzz_xxxy[k] = -g_y_0_yyyz_xxxy[k] * cd_z[k] + g_y_0_yyyz_xxxyz[k];

                g_y_0_yyyzz_xxxz[k] = -g_y_0_yyyz_xxxz[k] * cd_z[k] + g_y_0_yyyz_xxxzz[k];

                g_y_0_yyyzz_xxyy[k] = -g_y_0_yyyz_xxyy[k] * cd_z[k] + g_y_0_yyyz_xxyyz[k];

                g_y_0_yyyzz_xxyz[k] = -g_y_0_yyyz_xxyz[k] * cd_z[k] + g_y_0_yyyz_xxyzz[k];

                g_y_0_yyyzz_xxzz[k] = -g_y_0_yyyz_xxzz[k] * cd_z[k] + g_y_0_yyyz_xxzzz[k];

                g_y_0_yyyzz_xyyy[k] = -g_y_0_yyyz_xyyy[k] * cd_z[k] + g_y_0_yyyz_xyyyz[k];

                g_y_0_yyyzz_xyyz[k] = -g_y_0_yyyz_xyyz[k] * cd_z[k] + g_y_0_yyyz_xyyzz[k];

                g_y_0_yyyzz_xyzz[k] = -g_y_0_yyyz_xyzz[k] * cd_z[k] + g_y_0_yyyz_xyzzz[k];

                g_y_0_yyyzz_xzzz[k] = -g_y_0_yyyz_xzzz[k] * cd_z[k] + g_y_0_yyyz_xzzzz[k];

                g_y_0_yyyzz_yyyy[k] = -g_y_0_yyyz_yyyy[k] * cd_z[k] + g_y_0_yyyz_yyyyz[k];

                g_y_0_yyyzz_yyyz[k] = -g_y_0_yyyz_yyyz[k] * cd_z[k] + g_y_0_yyyz_yyyzz[k];

                g_y_0_yyyzz_yyzz[k] = -g_y_0_yyyz_yyzz[k] * cd_z[k] + g_y_0_yyyz_yyzzz[k];

                g_y_0_yyyzz_yzzz[k] = -g_y_0_yyyz_yzzz[k] * cd_z[k] + g_y_0_yyyz_yzzzz[k];

                g_y_0_yyyzz_zzzz[k] = -g_y_0_yyyz_zzzz[k] * cd_z[k] + g_y_0_yyyz_zzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 270);

            auto g_y_0_yyzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 271);

            auto g_y_0_yyzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 272);

            auto g_y_0_yyzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 273);

            auto g_y_0_yyzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 274);

            auto g_y_0_yyzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 275);

            auto g_y_0_yyzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 276);

            auto g_y_0_yyzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 277);

            auto g_y_0_yyzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 278);

            auto g_y_0_yyzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 279);

            auto g_y_0_yyzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 280);

            auto g_y_0_yyzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 281);

            auto g_y_0_yyzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 282);

            auto g_y_0_yyzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 283);

            auto g_y_0_yyzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 284);

            #pragma omp simd aligned(cd_z, g_y_0_yyzz_xxxx, g_y_0_yyzz_xxxxz, g_y_0_yyzz_xxxy, g_y_0_yyzz_xxxyz, g_y_0_yyzz_xxxz, g_y_0_yyzz_xxxzz, g_y_0_yyzz_xxyy, g_y_0_yyzz_xxyyz, g_y_0_yyzz_xxyz, g_y_0_yyzz_xxyzz, g_y_0_yyzz_xxzz, g_y_0_yyzz_xxzzz, g_y_0_yyzz_xyyy, g_y_0_yyzz_xyyyz, g_y_0_yyzz_xyyz, g_y_0_yyzz_xyyzz, g_y_0_yyzz_xyzz, g_y_0_yyzz_xyzzz, g_y_0_yyzz_xzzz, g_y_0_yyzz_xzzzz, g_y_0_yyzz_yyyy, g_y_0_yyzz_yyyyz, g_y_0_yyzz_yyyz, g_y_0_yyzz_yyyzz, g_y_0_yyzz_yyzz, g_y_0_yyzz_yyzzz, g_y_0_yyzz_yzzz, g_y_0_yyzz_yzzzz, g_y_0_yyzz_zzzz, g_y_0_yyzz_zzzzz, g_y_0_yyzzz_xxxx, g_y_0_yyzzz_xxxy, g_y_0_yyzzz_xxxz, g_y_0_yyzzz_xxyy, g_y_0_yyzzz_xxyz, g_y_0_yyzzz_xxzz, g_y_0_yyzzz_xyyy, g_y_0_yyzzz_xyyz, g_y_0_yyzzz_xyzz, g_y_0_yyzzz_xzzz, g_y_0_yyzzz_yyyy, g_y_0_yyzzz_yyyz, g_y_0_yyzzz_yyzz, g_y_0_yyzzz_yzzz, g_y_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_xxxx[k] = -g_y_0_yyzz_xxxx[k] * cd_z[k] + g_y_0_yyzz_xxxxz[k];

                g_y_0_yyzzz_xxxy[k] = -g_y_0_yyzz_xxxy[k] * cd_z[k] + g_y_0_yyzz_xxxyz[k];

                g_y_0_yyzzz_xxxz[k] = -g_y_0_yyzz_xxxz[k] * cd_z[k] + g_y_0_yyzz_xxxzz[k];

                g_y_0_yyzzz_xxyy[k] = -g_y_0_yyzz_xxyy[k] * cd_z[k] + g_y_0_yyzz_xxyyz[k];

                g_y_0_yyzzz_xxyz[k] = -g_y_0_yyzz_xxyz[k] * cd_z[k] + g_y_0_yyzz_xxyzz[k];

                g_y_0_yyzzz_xxzz[k] = -g_y_0_yyzz_xxzz[k] * cd_z[k] + g_y_0_yyzz_xxzzz[k];

                g_y_0_yyzzz_xyyy[k] = -g_y_0_yyzz_xyyy[k] * cd_z[k] + g_y_0_yyzz_xyyyz[k];

                g_y_0_yyzzz_xyyz[k] = -g_y_0_yyzz_xyyz[k] * cd_z[k] + g_y_0_yyzz_xyyzz[k];

                g_y_0_yyzzz_xyzz[k] = -g_y_0_yyzz_xyzz[k] * cd_z[k] + g_y_0_yyzz_xyzzz[k];

                g_y_0_yyzzz_xzzz[k] = -g_y_0_yyzz_xzzz[k] * cd_z[k] + g_y_0_yyzz_xzzzz[k];

                g_y_0_yyzzz_yyyy[k] = -g_y_0_yyzz_yyyy[k] * cd_z[k] + g_y_0_yyzz_yyyyz[k];

                g_y_0_yyzzz_yyyz[k] = -g_y_0_yyzz_yyyz[k] * cd_z[k] + g_y_0_yyzz_yyyzz[k];

                g_y_0_yyzzz_yyzz[k] = -g_y_0_yyzz_yyzz[k] * cd_z[k] + g_y_0_yyzz_yyzzz[k];

                g_y_0_yyzzz_yzzz[k] = -g_y_0_yyzz_yzzz[k] * cd_z[k] + g_y_0_yyzz_yzzzz[k];

                g_y_0_yyzzz_zzzz[k] = -g_y_0_yyzz_zzzz[k] * cd_z[k] + g_y_0_yyzz_zzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 285);

            auto g_y_0_yzzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 286);

            auto g_y_0_yzzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 287);

            auto g_y_0_yzzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 288);

            auto g_y_0_yzzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 289);

            auto g_y_0_yzzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 290);

            auto g_y_0_yzzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 291);

            auto g_y_0_yzzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 292);

            auto g_y_0_yzzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 293);

            auto g_y_0_yzzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 294);

            auto g_y_0_yzzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 295);

            auto g_y_0_yzzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 296);

            auto g_y_0_yzzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 297);

            auto g_y_0_yzzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 298);

            auto g_y_0_yzzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 299);

            #pragma omp simd aligned(cd_z, g_y_0_yzzz_xxxx, g_y_0_yzzz_xxxxz, g_y_0_yzzz_xxxy, g_y_0_yzzz_xxxyz, g_y_0_yzzz_xxxz, g_y_0_yzzz_xxxzz, g_y_0_yzzz_xxyy, g_y_0_yzzz_xxyyz, g_y_0_yzzz_xxyz, g_y_0_yzzz_xxyzz, g_y_0_yzzz_xxzz, g_y_0_yzzz_xxzzz, g_y_0_yzzz_xyyy, g_y_0_yzzz_xyyyz, g_y_0_yzzz_xyyz, g_y_0_yzzz_xyyzz, g_y_0_yzzz_xyzz, g_y_0_yzzz_xyzzz, g_y_0_yzzz_xzzz, g_y_0_yzzz_xzzzz, g_y_0_yzzz_yyyy, g_y_0_yzzz_yyyyz, g_y_0_yzzz_yyyz, g_y_0_yzzz_yyyzz, g_y_0_yzzz_yyzz, g_y_0_yzzz_yyzzz, g_y_0_yzzz_yzzz, g_y_0_yzzz_yzzzz, g_y_0_yzzz_zzzz, g_y_0_yzzz_zzzzz, g_y_0_yzzzz_xxxx, g_y_0_yzzzz_xxxy, g_y_0_yzzzz_xxxz, g_y_0_yzzzz_xxyy, g_y_0_yzzzz_xxyz, g_y_0_yzzzz_xxzz, g_y_0_yzzzz_xyyy, g_y_0_yzzzz_xyyz, g_y_0_yzzzz_xyzz, g_y_0_yzzzz_xzzz, g_y_0_yzzzz_yyyy, g_y_0_yzzzz_yyyz, g_y_0_yzzzz_yyzz, g_y_0_yzzzz_yzzz, g_y_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_xxxx[k] = -g_y_0_yzzz_xxxx[k] * cd_z[k] + g_y_0_yzzz_xxxxz[k];

                g_y_0_yzzzz_xxxy[k] = -g_y_0_yzzz_xxxy[k] * cd_z[k] + g_y_0_yzzz_xxxyz[k];

                g_y_0_yzzzz_xxxz[k] = -g_y_0_yzzz_xxxz[k] * cd_z[k] + g_y_0_yzzz_xxxzz[k];

                g_y_0_yzzzz_xxyy[k] = -g_y_0_yzzz_xxyy[k] * cd_z[k] + g_y_0_yzzz_xxyyz[k];

                g_y_0_yzzzz_xxyz[k] = -g_y_0_yzzz_xxyz[k] * cd_z[k] + g_y_0_yzzz_xxyzz[k];

                g_y_0_yzzzz_xxzz[k] = -g_y_0_yzzz_xxzz[k] * cd_z[k] + g_y_0_yzzz_xxzzz[k];

                g_y_0_yzzzz_xyyy[k] = -g_y_0_yzzz_xyyy[k] * cd_z[k] + g_y_0_yzzz_xyyyz[k];

                g_y_0_yzzzz_xyyz[k] = -g_y_0_yzzz_xyyz[k] * cd_z[k] + g_y_0_yzzz_xyyzz[k];

                g_y_0_yzzzz_xyzz[k] = -g_y_0_yzzz_xyzz[k] * cd_z[k] + g_y_0_yzzz_xyzzz[k];

                g_y_0_yzzzz_xzzz[k] = -g_y_0_yzzz_xzzz[k] * cd_z[k] + g_y_0_yzzz_xzzzz[k];

                g_y_0_yzzzz_yyyy[k] = -g_y_0_yzzz_yyyy[k] * cd_z[k] + g_y_0_yzzz_yyyyz[k];

                g_y_0_yzzzz_yyyz[k] = -g_y_0_yzzz_yyyz[k] * cd_z[k] + g_y_0_yzzz_yyyzz[k];

                g_y_0_yzzzz_yyzz[k] = -g_y_0_yzzz_yyzz[k] * cd_z[k] + g_y_0_yzzz_yyzzz[k];

                g_y_0_yzzzz_yzzz[k] = -g_y_0_yzzz_yzzz[k] * cd_z[k] + g_y_0_yzzz_yzzzz[k];

                g_y_0_yzzzz_zzzz[k] = -g_y_0_yzzz_zzzz[k] * cd_z[k] + g_y_0_yzzz_zzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 300);

            auto g_y_0_zzzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 301);

            auto g_y_0_zzzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 302);

            auto g_y_0_zzzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 303);

            auto g_y_0_zzzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 304);

            auto g_y_0_zzzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 305);

            auto g_y_0_zzzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 306);

            auto g_y_0_zzzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 307);

            auto g_y_0_zzzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 308);

            auto g_y_0_zzzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 309);

            auto g_y_0_zzzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 310);

            auto g_y_0_zzzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 311);

            auto g_y_0_zzzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 312);

            auto g_y_0_zzzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 313);

            auto g_y_0_zzzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_z, g_y_0_zzzz_xxxx, g_y_0_zzzz_xxxxz, g_y_0_zzzz_xxxy, g_y_0_zzzz_xxxyz, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxxzz, g_y_0_zzzz_xxyy, g_y_0_zzzz_xxyyz, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxyzz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xxzzz, g_y_0_zzzz_xyyy, g_y_0_zzzz_xyyyz, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyyzz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xyzzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_xzzzz, g_y_0_zzzz_yyyy, g_y_0_zzzz_yyyyz, g_y_0_zzzz_yyyz, g_y_0_zzzz_yyyzz, g_y_0_zzzz_yyzz, g_y_0_zzzz_yyzzz, g_y_0_zzzz_yzzz, g_y_0_zzzz_yzzzz, g_y_0_zzzz_zzzz, g_y_0_zzzz_zzzzz, g_y_0_zzzzz_xxxx, g_y_0_zzzzz_xxxy, g_y_0_zzzzz_xxxz, g_y_0_zzzzz_xxyy, g_y_0_zzzzz_xxyz, g_y_0_zzzzz_xxzz, g_y_0_zzzzz_xyyy, g_y_0_zzzzz_xyyz, g_y_0_zzzzz_xyzz, g_y_0_zzzzz_xzzz, g_y_0_zzzzz_yyyy, g_y_0_zzzzz_yyyz, g_y_0_zzzzz_yyzz, g_y_0_zzzzz_yzzz, g_y_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_xxxx[k] = -g_y_0_zzzz_xxxx[k] * cd_z[k] + g_y_0_zzzz_xxxxz[k];

                g_y_0_zzzzz_xxxy[k] = -g_y_0_zzzz_xxxy[k] * cd_z[k] + g_y_0_zzzz_xxxyz[k];

                g_y_0_zzzzz_xxxz[k] = -g_y_0_zzzz_xxxz[k] * cd_z[k] + g_y_0_zzzz_xxxzz[k];

                g_y_0_zzzzz_xxyy[k] = -g_y_0_zzzz_xxyy[k] * cd_z[k] + g_y_0_zzzz_xxyyz[k];

                g_y_0_zzzzz_xxyz[k] = -g_y_0_zzzz_xxyz[k] * cd_z[k] + g_y_0_zzzz_xxyzz[k];

                g_y_0_zzzzz_xxzz[k] = -g_y_0_zzzz_xxzz[k] * cd_z[k] + g_y_0_zzzz_xxzzz[k];

                g_y_0_zzzzz_xyyy[k] = -g_y_0_zzzz_xyyy[k] * cd_z[k] + g_y_0_zzzz_xyyyz[k];

                g_y_0_zzzzz_xyyz[k] = -g_y_0_zzzz_xyyz[k] * cd_z[k] + g_y_0_zzzz_xyyzz[k];

                g_y_0_zzzzz_xyzz[k] = -g_y_0_zzzz_xyzz[k] * cd_z[k] + g_y_0_zzzz_xyzzz[k];

                g_y_0_zzzzz_xzzz[k] = -g_y_0_zzzz_xzzz[k] * cd_z[k] + g_y_0_zzzz_xzzzz[k];

                g_y_0_zzzzz_yyyy[k] = -g_y_0_zzzz_yyyy[k] * cd_z[k] + g_y_0_zzzz_yyyyz[k];

                g_y_0_zzzzz_yyyz[k] = -g_y_0_zzzz_yyyz[k] * cd_z[k] + g_y_0_zzzz_yyyzz[k];

                g_y_0_zzzzz_yyzz[k] = -g_y_0_zzzz_yyzz[k] * cd_z[k] + g_y_0_zzzz_yyzzz[k];

                g_y_0_zzzzz_yzzz[k] = -g_y_0_zzzz_yzzz[k] * cd_z[k] + g_y_0_zzzz_yzzzz[k];

                g_y_0_zzzzz_zzzz[k] = -g_y_0_zzzz_zzzz[k] * cd_z[k] + g_y_0_zzzz_zzzzz[k];
            }
            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxx_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 5);

            auto g_z_0_xxxxx_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 6);

            auto g_z_0_xxxxx_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 7);

            auto g_z_0_xxxxx_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 8);

            auto g_z_0_xxxxx_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 9);

            auto g_z_0_xxxxx_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 10);

            auto g_z_0_xxxxx_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 11);

            auto g_z_0_xxxxx_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 12);

            auto g_z_0_xxxxx_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 13);

            auto g_z_0_xxxxx_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_z_0_xxxx_xxxx, g_z_0_xxxx_xxxxx, g_z_0_xxxx_xxxxy, g_z_0_xxxx_xxxxz, g_z_0_xxxx_xxxy, g_z_0_xxxx_xxxyy, g_z_0_xxxx_xxxyz, g_z_0_xxxx_xxxz, g_z_0_xxxx_xxxzz, g_z_0_xxxx_xxyy, g_z_0_xxxx_xxyyy, g_z_0_xxxx_xxyyz, g_z_0_xxxx_xxyz, g_z_0_xxxx_xxyzz, g_z_0_xxxx_xxzz, g_z_0_xxxx_xxzzz, g_z_0_xxxx_xyyy, g_z_0_xxxx_xyyyy, g_z_0_xxxx_xyyyz, g_z_0_xxxx_xyyz, g_z_0_xxxx_xyyzz, g_z_0_xxxx_xyzz, g_z_0_xxxx_xyzzz, g_z_0_xxxx_xzzz, g_z_0_xxxx_xzzzz, g_z_0_xxxx_yyyy, g_z_0_xxxx_yyyz, g_z_0_xxxx_yyzz, g_z_0_xxxx_yzzz, g_z_0_xxxx_zzzz, g_z_0_xxxxx_xxxx, g_z_0_xxxxx_xxxy, g_z_0_xxxxx_xxxz, g_z_0_xxxxx_xxyy, g_z_0_xxxxx_xxyz, g_z_0_xxxxx_xxzz, g_z_0_xxxxx_xyyy, g_z_0_xxxxx_xyyz, g_z_0_xxxxx_xyzz, g_z_0_xxxxx_xzzz, g_z_0_xxxxx_yyyy, g_z_0_xxxxx_yyyz, g_z_0_xxxxx_yyzz, g_z_0_xxxxx_yzzz, g_z_0_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_xxxx[k] = -g_z_0_xxxx_xxxx[k] * cd_x[k] + g_z_0_xxxx_xxxxx[k];

                g_z_0_xxxxx_xxxy[k] = -g_z_0_xxxx_xxxy[k] * cd_x[k] + g_z_0_xxxx_xxxxy[k];

                g_z_0_xxxxx_xxxz[k] = -g_z_0_xxxx_xxxz[k] * cd_x[k] + g_z_0_xxxx_xxxxz[k];

                g_z_0_xxxxx_xxyy[k] = -g_z_0_xxxx_xxyy[k] * cd_x[k] + g_z_0_xxxx_xxxyy[k];

                g_z_0_xxxxx_xxyz[k] = -g_z_0_xxxx_xxyz[k] * cd_x[k] + g_z_0_xxxx_xxxyz[k];

                g_z_0_xxxxx_xxzz[k] = -g_z_0_xxxx_xxzz[k] * cd_x[k] + g_z_0_xxxx_xxxzz[k];

                g_z_0_xxxxx_xyyy[k] = -g_z_0_xxxx_xyyy[k] * cd_x[k] + g_z_0_xxxx_xxyyy[k];

                g_z_0_xxxxx_xyyz[k] = -g_z_0_xxxx_xyyz[k] * cd_x[k] + g_z_0_xxxx_xxyyz[k];

                g_z_0_xxxxx_xyzz[k] = -g_z_0_xxxx_xyzz[k] * cd_x[k] + g_z_0_xxxx_xxyzz[k];

                g_z_0_xxxxx_xzzz[k] = -g_z_0_xxxx_xzzz[k] * cd_x[k] + g_z_0_xxxx_xxzzz[k];

                g_z_0_xxxxx_yyyy[k] = -g_z_0_xxxx_yyyy[k] * cd_x[k] + g_z_0_xxxx_xyyyy[k];

                g_z_0_xxxxx_yyyz[k] = -g_z_0_xxxx_yyyz[k] * cd_x[k] + g_z_0_xxxx_xyyyz[k];

                g_z_0_xxxxx_yyzz[k] = -g_z_0_xxxx_yyzz[k] * cd_x[k] + g_z_0_xxxx_xyyzz[k];

                g_z_0_xxxxx_yzzz[k] = -g_z_0_xxxx_yzzz[k] * cd_x[k] + g_z_0_xxxx_xyzzz[k];

                g_z_0_xxxxx_zzzz[k] = -g_z_0_xxxx_zzzz[k] * cd_x[k] + g_z_0_xxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxy_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 15);

            auto g_z_0_xxxxy_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 16);

            auto g_z_0_xxxxy_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 17);

            auto g_z_0_xxxxy_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 18);

            auto g_z_0_xxxxy_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 19);

            auto g_z_0_xxxxy_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 20);

            auto g_z_0_xxxxy_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 21);

            auto g_z_0_xxxxy_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 22);

            auto g_z_0_xxxxy_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 23);

            auto g_z_0_xxxxy_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 24);

            auto g_z_0_xxxxy_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 25);

            auto g_z_0_xxxxy_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 26);

            auto g_z_0_xxxxy_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 27);

            auto g_z_0_xxxxy_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 28);

            auto g_z_0_xxxxy_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxy_xxxx, g_z_0_xxxxy_xxxy, g_z_0_xxxxy_xxxz, g_z_0_xxxxy_xxyy, g_z_0_xxxxy_xxyz, g_z_0_xxxxy_xxzz, g_z_0_xxxxy_xyyy, g_z_0_xxxxy_xyyz, g_z_0_xxxxy_xyzz, g_z_0_xxxxy_xzzz, g_z_0_xxxxy_yyyy, g_z_0_xxxxy_yyyz, g_z_0_xxxxy_yyzz, g_z_0_xxxxy_yzzz, g_z_0_xxxxy_zzzz, g_z_0_xxxy_xxxx, g_z_0_xxxy_xxxxx, g_z_0_xxxy_xxxxy, g_z_0_xxxy_xxxxz, g_z_0_xxxy_xxxy, g_z_0_xxxy_xxxyy, g_z_0_xxxy_xxxyz, g_z_0_xxxy_xxxz, g_z_0_xxxy_xxxzz, g_z_0_xxxy_xxyy, g_z_0_xxxy_xxyyy, g_z_0_xxxy_xxyyz, g_z_0_xxxy_xxyz, g_z_0_xxxy_xxyzz, g_z_0_xxxy_xxzz, g_z_0_xxxy_xxzzz, g_z_0_xxxy_xyyy, g_z_0_xxxy_xyyyy, g_z_0_xxxy_xyyyz, g_z_0_xxxy_xyyz, g_z_0_xxxy_xyyzz, g_z_0_xxxy_xyzz, g_z_0_xxxy_xyzzz, g_z_0_xxxy_xzzz, g_z_0_xxxy_xzzzz, g_z_0_xxxy_yyyy, g_z_0_xxxy_yyyz, g_z_0_xxxy_yyzz, g_z_0_xxxy_yzzz, g_z_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_xxxx[k] = -g_z_0_xxxy_xxxx[k] * cd_x[k] + g_z_0_xxxy_xxxxx[k];

                g_z_0_xxxxy_xxxy[k] = -g_z_0_xxxy_xxxy[k] * cd_x[k] + g_z_0_xxxy_xxxxy[k];

                g_z_0_xxxxy_xxxz[k] = -g_z_0_xxxy_xxxz[k] * cd_x[k] + g_z_0_xxxy_xxxxz[k];

                g_z_0_xxxxy_xxyy[k] = -g_z_0_xxxy_xxyy[k] * cd_x[k] + g_z_0_xxxy_xxxyy[k];

                g_z_0_xxxxy_xxyz[k] = -g_z_0_xxxy_xxyz[k] * cd_x[k] + g_z_0_xxxy_xxxyz[k];

                g_z_0_xxxxy_xxzz[k] = -g_z_0_xxxy_xxzz[k] * cd_x[k] + g_z_0_xxxy_xxxzz[k];

                g_z_0_xxxxy_xyyy[k] = -g_z_0_xxxy_xyyy[k] * cd_x[k] + g_z_0_xxxy_xxyyy[k];

                g_z_0_xxxxy_xyyz[k] = -g_z_0_xxxy_xyyz[k] * cd_x[k] + g_z_0_xxxy_xxyyz[k];

                g_z_0_xxxxy_xyzz[k] = -g_z_0_xxxy_xyzz[k] * cd_x[k] + g_z_0_xxxy_xxyzz[k];

                g_z_0_xxxxy_xzzz[k] = -g_z_0_xxxy_xzzz[k] * cd_x[k] + g_z_0_xxxy_xxzzz[k];

                g_z_0_xxxxy_yyyy[k] = -g_z_0_xxxy_yyyy[k] * cd_x[k] + g_z_0_xxxy_xyyyy[k];

                g_z_0_xxxxy_yyyz[k] = -g_z_0_xxxy_yyyz[k] * cd_x[k] + g_z_0_xxxy_xyyyz[k];

                g_z_0_xxxxy_yyzz[k] = -g_z_0_xxxy_yyzz[k] * cd_x[k] + g_z_0_xxxy_xyyzz[k];

                g_z_0_xxxxy_yzzz[k] = -g_z_0_xxxy_yzzz[k] * cd_x[k] + g_z_0_xxxy_xyzzz[k];

                g_z_0_xxxxy_zzzz[k] = -g_z_0_xxxy_zzzz[k] * cd_x[k] + g_z_0_xxxy_xzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 30);

            auto g_z_0_xxxxz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 31);

            auto g_z_0_xxxxz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 32);

            auto g_z_0_xxxxz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 33);

            auto g_z_0_xxxxz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 34);

            auto g_z_0_xxxxz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 35);

            auto g_z_0_xxxxz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 36);

            auto g_z_0_xxxxz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 37);

            auto g_z_0_xxxxz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 38);

            auto g_z_0_xxxxz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 39);

            auto g_z_0_xxxxz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 40);

            auto g_z_0_xxxxz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 41);

            auto g_z_0_xxxxz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 42);

            auto g_z_0_xxxxz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 43);

            auto g_z_0_xxxxz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxz_xxxx, g_z_0_xxxxz_xxxy, g_z_0_xxxxz_xxxz, g_z_0_xxxxz_xxyy, g_z_0_xxxxz_xxyz, g_z_0_xxxxz_xxzz, g_z_0_xxxxz_xyyy, g_z_0_xxxxz_xyyz, g_z_0_xxxxz_xyzz, g_z_0_xxxxz_xzzz, g_z_0_xxxxz_yyyy, g_z_0_xxxxz_yyyz, g_z_0_xxxxz_yyzz, g_z_0_xxxxz_yzzz, g_z_0_xxxxz_zzzz, g_z_0_xxxz_xxxx, g_z_0_xxxz_xxxxx, g_z_0_xxxz_xxxxy, g_z_0_xxxz_xxxxz, g_z_0_xxxz_xxxy, g_z_0_xxxz_xxxyy, g_z_0_xxxz_xxxyz, g_z_0_xxxz_xxxz, g_z_0_xxxz_xxxzz, g_z_0_xxxz_xxyy, g_z_0_xxxz_xxyyy, g_z_0_xxxz_xxyyz, g_z_0_xxxz_xxyz, g_z_0_xxxz_xxyzz, g_z_0_xxxz_xxzz, g_z_0_xxxz_xxzzz, g_z_0_xxxz_xyyy, g_z_0_xxxz_xyyyy, g_z_0_xxxz_xyyyz, g_z_0_xxxz_xyyz, g_z_0_xxxz_xyyzz, g_z_0_xxxz_xyzz, g_z_0_xxxz_xyzzz, g_z_0_xxxz_xzzz, g_z_0_xxxz_xzzzz, g_z_0_xxxz_yyyy, g_z_0_xxxz_yyyz, g_z_0_xxxz_yyzz, g_z_0_xxxz_yzzz, g_z_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_xxxx[k] = -g_z_0_xxxz_xxxx[k] * cd_x[k] + g_z_0_xxxz_xxxxx[k];

                g_z_0_xxxxz_xxxy[k] = -g_z_0_xxxz_xxxy[k] * cd_x[k] + g_z_0_xxxz_xxxxy[k];

                g_z_0_xxxxz_xxxz[k] = -g_z_0_xxxz_xxxz[k] * cd_x[k] + g_z_0_xxxz_xxxxz[k];

                g_z_0_xxxxz_xxyy[k] = -g_z_0_xxxz_xxyy[k] * cd_x[k] + g_z_0_xxxz_xxxyy[k];

                g_z_0_xxxxz_xxyz[k] = -g_z_0_xxxz_xxyz[k] * cd_x[k] + g_z_0_xxxz_xxxyz[k];

                g_z_0_xxxxz_xxzz[k] = -g_z_0_xxxz_xxzz[k] * cd_x[k] + g_z_0_xxxz_xxxzz[k];

                g_z_0_xxxxz_xyyy[k] = -g_z_0_xxxz_xyyy[k] * cd_x[k] + g_z_0_xxxz_xxyyy[k];

                g_z_0_xxxxz_xyyz[k] = -g_z_0_xxxz_xyyz[k] * cd_x[k] + g_z_0_xxxz_xxyyz[k];

                g_z_0_xxxxz_xyzz[k] = -g_z_0_xxxz_xyzz[k] * cd_x[k] + g_z_0_xxxz_xxyzz[k];

                g_z_0_xxxxz_xzzz[k] = -g_z_0_xxxz_xzzz[k] * cd_x[k] + g_z_0_xxxz_xxzzz[k];

                g_z_0_xxxxz_yyyy[k] = -g_z_0_xxxz_yyyy[k] * cd_x[k] + g_z_0_xxxz_xyyyy[k];

                g_z_0_xxxxz_yyyz[k] = -g_z_0_xxxz_yyyz[k] * cd_x[k] + g_z_0_xxxz_xyyyz[k];

                g_z_0_xxxxz_yyzz[k] = -g_z_0_xxxz_yyzz[k] * cd_x[k] + g_z_0_xxxz_xyyzz[k];

                g_z_0_xxxxz_yzzz[k] = -g_z_0_xxxz_yzzz[k] * cd_x[k] + g_z_0_xxxz_xyzzz[k];

                g_z_0_xxxxz_zzzz[k] = -g_z_0_xxxz_zzzz[k] * cd_x[k] + g_z_0_xxxz_xzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyy_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 45);

            auto g_z_0_xxxyy_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 46);

            auto g_z_0_xxxyy_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 47);

            auto g_z_0_xxxyy_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 48);

            auto g_z_0_xxxyy_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 49);

            auto g_z_0_xxxyy_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 50);

            auto g_z_0_xxxyy_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 51);

            auto g_z_0_xxxyy_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 52);

            auto g_z_0_xxxyy_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 53);

            auto g_z_0_xxxyy_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 54);

            auto g_z_0_xxxyy_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 55);

            auto g_z_0_xxxyy_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 56);

            auto g_z_0_xxxyy_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 57);

            auto g_z_0_xxxyy_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 58);

            auto g_z_0_xxxyy_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyy_xxxx, g_z_0_xxxyy_xxxy, g_z_0_xxxyy_xxxz, g_z_0_xxxyy_xxyy, g_z_0_xxxyy_xxyz, g_z_0_xxxyy_xxzz, g_z_0_xxxyy_xyyy, g_z_0_xxxyy_xyyz, g_z_0_xxxyy_xyzz, g_z_0_xxxyy_xzzz, g_z_0_xxxyy_yyyy, g_z_0_xxxyy_yyyz, g_z_0_xxxyy_yyzz, g_z_0_xxxyy_yzzz, g_z_0_xxxyy_zzzz, g_z_0_xxyy_xxxx, g_z_0_xxyy_xxxxx, g_z_0_xxyy_xxxxy, g_z_0_xxyy_xxxxz, g_z_0_xxyy_xxxy, g_z_0_xxyy_xxxyy, g_z_0_xxyy_xxxyz, g_z_0_xxyy_xxxz, g_z_0_xxyy_xxxzz, g_z_0_xxyy_xxyy, g_z_0_xxyy_xxyyy, g_z_0_xxyy_xxyyz, g_z_0_xxyy_xxyz, g_z_0_xxyy_xxyzz, g_z_0_xxyy_xxzz, g_z_0_xxyy_xxzzz, g_z_0_xxyy_xyyy, g_z_0_xxyy_xyyyy, g_z_0_xxyy_xyyyz, g_z_0_xxyy_xyyz, g_z_0_xxyy_xyyzz, g_z_0_xxyy_xyzz, g_z_0_xxyy_xyzzz, g_z_0_xxyy_xzzz, g_z_0_xxyy_xzzzz, g_z_0_xxyy_yyyy, g_z_0_xxyy_yyyz, g_z_0_xxyy_yyzz, g_z_0_xxyy_yzzz, g_z_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_xxxx[k] = -g_z_0_xxyy_xxxx[k] * cd_x[k] + g_z_0_xxyy_xxxxx[k];

                g_z_0_xxxyy_xxxy[k] = -g_z_0_xxyy_xxxy[k] * cd_x[k] + g_z_0_xxyy_xxxxy[k];

                g_z_0_xxxyy_xxxz[k] = -g_z_0_xxyy_xxxz[k] * cd_x[k] + g_z_0_xxyy_xxxxz[k];

                g_z_0_xxxyy_xxyy[k] = -g_z_0_xxyy_xxyy[k] * cd_x[k] + g_z_0_xxyy_xxxyy[k];

                g_z_0_xxxyy_xxyz[k] = -g_z_0_xxyy_xxyz[k] * cd_x[k] + g_z_0_xxyy_xxxyz[k];

                g_z_0_xxxyy_xxzz[k] = -g_z_0_xxyy_xxzz[k] * cd_x[k] + g_z_0_xxyy_xxxzz[k];

                g_z_0_xxxyy_xyyy[k] = -g_z_0_xxyy_xyyy[k] * cd_x[k] + g_z_0_xxyy_xxyyy[k];

                g_z_0_xxxyy_xyyz[k] = -g_z_0_xxyy_xyyz[k] * cd_x[k] + g_z_0_xxyy_xxyyz[k];

                g_z_0_xxxyy_xyzz[k] = -g_z_0_xxyy_xyzz[k] * cd_x[k] + g_z_0_xxyy_xxyzz[k];

                g_z_0_xxxyy_xzzz[k] = -g_z_0_xxyy_xzzz[k] * cd_x[k] + g_z_0_xxyy_xxzzz[k];

                g_z_0_xxxyy_yyyy[k] = -g_z_0_xxyy_yyyy[k] * cd_x[k] + g_z_0_xxyy_xyyyy[k];

                g_z_0_xxxyy_yyyz[k] = -g_z_0_xxyy_yyyz[k] * cd_x[k] + g_z_0_xxyy_xyyyz[k];

                g_z_0_xxxyy_yyzz[k] = -g_z_0_xxyy_yyzz[k] * cd_x[k] + g_z_0_xxyy_xyyzz[k];

                g_z_0_xxxyy_yzzz[k] = -g_z_0_xxyy_yzzz[k] * cd_x[k] + g_z_0_xxyy_xyzzz[k];

                g_z_0_xxxyy_zzzz[k] = -g_z_0_xxyy_zzzz[k] * cd_x[k] + g_z_0_xxyy_xzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 60);

            auto g_z_0_xxxyz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 61);

            auto g_z_0_xxxyz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 62);

            auto g_z_0_xxxyz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 63);

            auto g_z_0_xxxyz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 64);

            auto g_z_0_xxxyz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 65);

            auto g_z_0_xxxyz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 66);

            auto g_z_0_xxxyz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 67);

            auto g_z_0_xxxyz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 68);

            auto g_z_0_xxxyz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 69);

            auto g_z_0_xxxyz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 70);

            auto g_z_0_xxxyz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 71);

            auto g_z_0_xxxyz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 72);

            auto g_z_0_xxxyz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 73);

            auto g_z_0_xxxyz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 74);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyz_xxxx, g_z_0_xxxyz_xxxy, g_z_0_xxxyz_xxxz, g_z_0_xxxyz_xxyy, g_z_0_xxxyz_xxyz, g_z_0_xxxyz_xxzz, g_z_0_xxxyz_xyyy, g_z_0_xxxyz_xyyz, g_z_0_xxxyz_xyzz, g_z_0_xxxyz_xzzz, g_z_0_xxxyz_yyyy, g_z_0_xxxyz_yyyz, g_z_0_xxxyz_yyzz, g_z_0_xxxyz_yzzz, g_z_0_xxxyz_zzzz, g_z_0_xxyz_xxxx, g_z_0_xxyz_xxxxx, g_z_0_xxyz_xxxxy, g_z_0_xxyz_xxxxz, g_z_0_xxyz_xxxy, g_z_0_xxyz_xxxyy, g_z_0_xxyz_xxxyz, g_z_0_xxyz_xxxz, g_z_0_xxyz_xxxzz, g_z_0_xxyz_xxyy, g_z_0_xxyz_xxyyy, g_z_0_xxyz_xxyyz, g_z_0_xxyz_xxyz, g_z_0_xxyz_xxyzz, g_z_0_xxyz_xxzz, g_z_0_xxyz_xxzzz, g_z_0_xxyz_xyyy, g_z_0_xxyz_xyyyy, g_z_0_xxyz_xyyyz, g_z_0_xxyz_xyyz, g_z_0_xxyz_xyyzz, g_z_0_xxyz_xyzz, g_z_0_xxyz_xyzzz, g_z_0_xxyz_xzzz, g_z_0_xxyz_xzzzz, g_z_0_xxyz_yyyy, g_z_0_xxyz_yyyz, g_z_0_xxyz_yyzz, g_z_0_xxyz_yzzz, g_z_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_xxxx[k] = -g_z_0_xxyz_xxxx[k] * cd_x[k] + g_z_0_xxyz_xxxxx[k];

                g_z_0_xxxyz_xxxy[k] = -g_z_0_xxyz_xxxy[k] * cd_x[k] + g_z_0_xxyz_xxxxy[k];

                g_z_0_xxxyz_xxxz[k] = -g_z_0_xxyz_xxxz[k] * cd_x[k] + g_z_0_xxyz_xxxxz[k];

                g_z_0_xxxyz_xxyy[k] = -g_z_0_xxyz_xxyy[k] * cd_x[k] + g_z_0_xxyz_xxxyy[k];

                g_z_0_xxxyz_xxyz[k] = -g_z_0_xxyz_xxyz[k] * cd_x[k] + g_z_0_xxyz_xxxyz[k];

                g_z_0_xxxyz_xxzz[k] = -g_z_0_xxyz_xxzz[k] * cd_x[k] + g_z_0_xxyz_xxxzz[k];

                g_z_0_xxxyz_xyyy[k] = -g_z_0_xxyz_xyyy[k] * cd_x[k] + g_z_0_xxyz_xxyyy[k];

                g_z_0_xxxyz_xyyz[k] = -g_z_0_xxyz_xyyz[k] * cd_x[k] + g_z_0_xxyz_xxyyz[k];

                g_z_0_xxxyz_xyzz[k] = -g_z_0_xxyz_xyzz[k] * cd_x[k] + g_z_0_xxyz_xxyzz[k];

                g_z_0_xxxyz_xzzz[k] = -g_z_0_xxyz_xzzz[k] * cd_x[k] + g_z_0_xxyz_xxzzz[k];

                g_z_0_xxxyz_yyyy[k] = -g_z_0_xxyz_yyyy[k] * cd_x[k] + g_z_0_xxyz_xyyyy[k];

                g_z_0_xxxyz_yyyz[k] = -g_z_0_xxyz_yyyz[k] * cd_x[k] + g_z_0_xxyz_xyyyz[k];

                g_z_0_xxxyz_yyzz[k] = -g_z_0_xxyz_yyzz[k] * cd_x[k] + g_z_0_xxyz_xyyzz[k];

                g_z_0_xxxyz_yzzz[k] = -g_z_0_xxyz_yzzz[k] * cd_x[k] + g_z_0_xxyz_xyzzz[k];

                g_z_0_xxxyz_zzzz[k] = -g_z_0_xxyz_zzzz[k] * cd_x[k] + g_z_0_xxyz_xzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 75);

            auto g_z_0_xxxzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 76);

            auto g_z_0_xxxzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 77);

            auto g_z_0_xxxzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 78);

            auto g_z_0_xxxzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 79);

            auto g_z_0_xxxzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 80);

            auto g_z_0_xxxzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 81);

            auto g_z_0_xxxzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 82);

            auto g_z_0_xxxzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 83);

            auto g_z_0_xxxzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 84);

            auto g_z_0_xxxzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 85);

            auto g_z_0_xxxzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 86);

            auto g_z_0_xxxzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 87);

            auto g_z_0_xxxzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 88);

            auto g_z_0_xxxzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzz_xxxx, g_z_0_xxxzz_xxxy, g_z_0_xxxzz_xxxz, g_z_0_xxxzz_xxyy, g_z_0_xxxzz_xxyz, g_z_0_xxxzz_xxzz, g_z_0_xxxzz_xyyy, g_z_0_xxxzz_xyyz, g_z_0_xxxzz_xyzz, g_z_0_xxxzz_xzzz, g_z_0_xxxzz_yyyy, g_z_0_xxxzz_yyyz, g_z_0_xxxzz_yyzz, g_z_0_xxxzz_yzzz, g_z_0_xxxzz_zzzz, g_z_0_xxzz_xxxx, g_z_0_xxzz_xxxxx, g_z_0_xxzz_xxxxy, g_z_0_xxzz_xxxxz, g_z_0_xxzz_xxxy, g_z_0_xxzz_xxxyy, g_z_0_xxzz_xxxyz, g_z_0_xxzz_xxxz, g_z_0_xxzz_xxxzz, g_z_0_xxzz_xxyy, g_z_0_xxzz_xxyyy, g_z_0_xxzz_xxyyz, g_z_0_xxzz_xxyz, g_z_0_xxzz_xxyzz, g_z_0_xxzz_xxzz, g_z_0_xxzz_xxzzz, g_z_0_xxzz_xyyy, g_z_0_xxzz_xyyyy, g_z_0_xxzz_xyyyz, g_z_0_xxzz_xyyz, g_z_0_xxzz_xyyzz, g_z_0_xxzz_xyzz, g_z_0_xxzz_xyzzz, g_z_0_xxzz_xzzz, g_z_0_xxzz_xzzzz, g_z_0_xxzz_yyyy, g_z_0_xxzz_yyyz, g_z_0_xxzz_yyzz, g_z_0_xxzz_yzzz, g_z_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_xxxx[k] = -g_z_0_xxzz_xxxx[k] * cd_x[k] + g_z_0_xxzz_xxxxx[k];

                g_z_0_xxxzz_xxxy[k] = -g_z_0_xxzz_xxxy[k] * cd_x[k] + g_z_0_xxzz_xxxxy[k];

                g_z_0_xxxzz_xxxz[k] = -g_z_0_xxzz_xxxz[k] * cd_x[k] + g_z_0_xxzz_xxxxz[k];

                g_z_0_xxxzz_xxyy[k] = -g_z_0_xxzz_xxyy[k] * cd_x[k] + g_z_0_xxzz_xxxyy[k];

                g_z_0_xxxzz_xxyz[k] = -g_z_0_xxzz_xxyz[k] * cd_x[k] + g_z_0_xxzz_xxxyz[k];

                g_z_0_xxxzz_xxzz[k] = -g_z_0_xxzz_xxzz[k] * cd_x[k] + g_z_0_xxzz_xxxzz[k];

                g_z_0_xxxzz_xyyy[k] = -g_z_0_xxzz_xyyy[k] * cd_x[k] + g_z_0_xxzz_xxyyy[k];

                g_z_0_xxxzz_xyyz[k] = -g_z_0_xxzz_xyyz[k] * cd_x[k] + g_z_0_xxzz_xxyyz[k];

                g_z_0_xxxzz_xyzz[k] = -g_z_0_xxzz_xyzz[k] * cd_x[k] + g_z_0_xxzz_xxyzz[k];

                g_z_0_xxxzz_xzzz[k] = -g_z_0_xxzz_xzzz[k] * cd_x[k] + g_z_0_xxzz_xxzzz[k];

                g_z_0_xxxzz_yyyy[k] = -g_z_0_xxzz_yyyy[k] * cd_x[k] + g_z_0_xxzz_xyyyy[k];

                g_z_0_xxxzz_yyyz[k] = -g_z_0_xxzz_yyyz[k] * cd_x[k] + g_z_0_xxzz_xyyyz[k];

                g_z_0_xxxzz_yyzz[k] = -g_z_0_xxzz_yyzz[k] * cd_x[k] + g_z_0_xxzz_xyyzz[k];

                g_z_0_xxxzz_yzzz[k] = -g_z_0_xxzz_yzzz[k] * cd_x[k] + g_z_0_xxzz_xyzzz[k];

                g_z_0_xxxzz_zzzz[k] = -g_z_0_xxzz_zzzz[k] * cd_x[k] + g_z_0_xxzz_xzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyy_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 90);

            auto g_z_0_xxyyy_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 91);

            auto g_z_0_xxyyy_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 92);

            auto g_z_0_xxyyy_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 93);

            auto g_z_0_xxyyy_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 94);

            auto g_z_0_xxyyy_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 95);

            auto g_z_0_xxyyy_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 96);

            auto g_z_0_xxyyy_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 97);

            auto g_z_0_xxyyy_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 98);

            auto g_z_0_xxyyy_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 99);

            auto g_z_0_xxyyy_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 100);

            auto g_z_0_xxyyy_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 101);

            auto g_z_0_xxyyy_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 102);

            auto g_z_0_xxyyy_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 103);

            auto g_z_0_xxyyy_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyy_xxxx, g_z_0_xxyyy_xxxy, g_z_0_xxyyy_xxxz, g_z_0_xxyyy_xxyy, g_z_0_xxyyy_xxyz, g_z_0_xxyyy_xxzz, g_z_0_xxyyy_xyyy, g_z_0_xxyyy_xyyz, g_z_0_xxyyy_xyzz, g_z_0_xxyyy_xzzz, g_z_0_xxyyy_yyyy, g_z_0_xxyyy_yyyz, g_z_0_xxyyy_yyzz, g_z_0_xxyyy_yzzz, g_z_0_xxyyy_zzzz, g_z_0_xyyy_xxxx, g_z_0_xyyy_xxxxx, g_z_0_xyyy_xxxxy, g_z_0_xyyy_xxxxz, g_z_0_xyyy_xxxy, g_z_0_xyyy_xxxyy, g_z_0_xyyy_xxxyz, g_z_0_xyyy_xxxz, g_z_0_xyyy_xxxzz, g_z_0_xyyy_xxyy, g_z_0_xyyy_xxyyy, g_z_0_xyyy_xxyyz, g_z_0_xyyy_xxyz, g_z_0_xyyy_xxyzz, g_z_0_xyyy_xxzz, g_z_0_xyyy_xxzzz, g_z_0_xyyy_xyyy, g_z_0_xyyy_xyyyy, g_z_0_xyyy_xyyyz, g_z_0_xyyy_xyyz, g_z_0_xyyy_xyyzz, g_z_0_xyyy_xyzz, g_z_0_xyyy_xyzzz, g_z_0_xyyy_xzzz, g_z_0_xyyy_xzzzz, g_z_0_xyyy_yyyy, g_z_0_xyyy_yyyz, g_z_0_xyyy_yyzz, g_z_0_xyyy_yzzz, g_z_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_xxxx[k] = -g_z_0_xyyy_xxxx[k] * cd_x[k] + g_z_0_xyyy_xxxxx[k];

                g_z_0_xxyyy_xxxy[k] = -g_z_0_xyyy_xxxy[k] * cd_x[k] + g_z_0_xyyy_xxxxy[k];

                g_z_0_xxyyy_xxxz[k] = -g_z_0_xyyy_xxxz[k] * cd_x[k] + g_z_0_xyyy_xxxxz[k];

                g_z_0_xxyyy_xxyy[k] = -g_z_0_xyyy_xxyy[k] * cd_x[k] + g_z_0_xyyy_xxxyy[k];

                g_z_0_xxyyy_xxyz[k] = -g_z_0_xyyy_xxyz[k] * cd_x[k] + g_z_0_xyyy_xxxyz[k];

                g_z_0_xxyyy_xxzz[k] = -g_z_0_xyyy_xxzz[k] * cd_x[k] + g_z_0_xyyy_xxxzz[k];

                g_z_0_xxyyy_xyyy[k] = -g_z_0_xyyy_xyyy[k] * cd_x[k] + g_z_0_xyyy_xxyyy[k];

                g_z_0_xxyyy_xyyz[k] = -g_z_0_xyyy_xyyz[k] * cd_x[k] + g_z_0_xyyy_xxyyz[k];

                g_z_0_xxyyy_xyzz[k] = -g_z_0_xyyy_xyzz[k] * cd_x[k] + g_z_0_xyyy_xxyzz[k];

                g_z_0_xxyyy_xzzz[k] = -g_z_0_xyyy_xzzz[k] * cd_x[k] + g_z_0_xyyy_xxzzz[k];

                g_z_0_xxyyy_yyyy[k] = -g_z_0_xyyy_yyyy[k] * cd_x[k] + g_z_0_xyyy_xyyyy[k];

                g_z_0_xxyyy_yyyz[k] = -g_z_0_xyyy_yyyz[k] * cd_x[k] + g_z_0_xyyy_xyyyz[k];

                g_z_0_xxyyy_yyzz[k] = -g_z_0_xyyy_yyzz[k] * cd_x[k] + g_z_0_xyyy_xyyzz[k];

                g_z_0_xxyyy_yzzz[k] = -g_z_0_xyyy_yzzz[k] * cd_x[k] + g_z_0_xyyy_xyzzz[k];

                g_z_0_xxyyy_zzzz[k] = -g_z_0_xyyy_zzzz[k] * cd_x[k] + g_z_0_xyyy_xzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 105);

            auto g_z_0_xxyyz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 106);

            auto g_z_0_xxyyz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 107);

            auto g_z_0_xxyyz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 108);

            auto g_z_0_xxyyz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 109);

            auto g_z_0_xxyyz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 110);

            auto g_z_0_xxyyz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 111);

            auto g_z_0_xxyyz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 112);

            auto g_z_0_xxyyz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 113);

            auto g_z_0_xxyyz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 114);

            auto g_z_0_xxyyz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 115);

            auto g_z_0_xxyyz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 116);

            auto g_z_0_xxyyz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 117);

            auto g_z_0_xxyyz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 118);

            auto g_z_0_xxyyz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyz_xxxx, g_z_0_xxyyz_xxxy, g_z_0_xxyyz_xxxz, g_z_0_xxyyz_xxyy, g_z_0_xxyyz_xxyz, g_z_0_xxyyz_xxzz, g_z_0_xxyyz_xyyy, g_z_0_xxyyz_xyyz, g_z_0_xxyyz_xyzz, g_z_0_xxyyz_xzzz, g_z_0_xxyyz_yyyy, g_z_0_xxyyz_yyyz, g_z_0_xxyyz_yyzz, g_z_0_xxyyz_yzzz, g_z_0_xxyyz_zzzz, g_z_0_xyyz_xxxx, g_z_0_xyyz_xxxxx, g_z_0_xyyz_xxxxy, g_z_0_xyyz_xxxxz, g_z_0_xyyz_xxxy, g_z_0_xyyz_xxxyy, g_z_0_xyyz_xxxyz, g_z_0_xyyz_xxxz, g_z_0_xyyz_xxxzz, g_z_0_xyyz_xxyy, g_z_0_xyyz_xxyyy, g_z_0_xyyz_xxyyz, g_z_0_xyyz_xxyz, g_z_0_xyyz_xxyzz, g_z_0_xyyz_xxzz, g_z_0_xyyz_xxzzz, g_z_0_xyyz_xyyy, g_z_0_xyyz_xyyyy, g_z_0_xyyz_xyyyz, g_z_0_xyyz_xyyz, g_z_0_xyyz_xyyzz, g_z_0_xyyz_xyzz, g_z_0_xyyz_xyzzz, g_z_0_xyyz_xzzz, g_z_0_xyyz_xzzzz, g_z_0_xyyz_yyyy, g_z_0_xyyz_yyyz, g_z_0_xyyz_yyzz, g_z_0_xyyz_yzzz, g_z_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_xxxx[k] = -g_z_0_xyyz_xxxx[k] * cd_x[k] + g_z_0_xyyz_xxxxx[k];

                g_z_0_xxyyz_xxxy[k] = -g_z_0_xyyz_xxxy[k] * cd_x[k] + g_z_0_xyyz_xxxxy[k];

                g_z_0_xxyyz_xxxz[k] = -g_z_0_xyyz_xxxz[k] * cd_x[k] + g_z_0_xyyz_xxxxz[k];

                g_z_0_xxyyz_xxyy[k] = -g_z_0_xyyz_xxyy[k] * cd_x[k] + g_z_0_xyyz_xxxyy[k];

                g_z_0_xxyyz_xxyz[k] = -g_z_0_xyyz_xxyz[k] * cd_x[k] + g_z_0_xyyz_xxxyz[k];

                g_z_0_xxyyz_xxzz[k] = -g_z_0_xyyz_xxzz[k] * cd_x[k] + g_z_0_xyyz_xxxzz[k];

                g_z_0_xxyyz_xyyy[k] = -g_z_0_xyyz_xyyy[k] * cd_x[k] + g_z_0_xyyz_xxyyy[k];

                g_z_0_xxyyz_xyyz[k] = -g_z_0_xyyz_xyyz[k] * cd_x[k] + g_z_0_xyyz_xxyyz[k];

                g_z_0_xxyyz_xyzz[k] = -g_z_0_xyyz_xyzz[k] * cd_x[k] + g_z_0_xyyz_xxyzz[k];

                g_z_0_xxyyz_xzzz[k] = -g_z_0_xyyz_xzzz[k] * cd_x[k] + g_z_0_xyyz_xxzzz[k];

                g_z_0_xxyyz_yyyy[k] = -g_z_0_xyyz_yyyy[k] * cd_x[k] + g_z_0_xyyz_xyyyy[k];

                g_z_0_xxyyz_yyyz[k] = -g_z_0_xyyz_yyyz[k] * cd_x[k] + g_z_0_xyyz_xyyyz[k];

                g_z_0_xxyyz_yyzz[k] = -g_z_0_xyyz_yyzz[k] * cd_x[k] + g_z_0_xyyz_xyyzz[k];

                g_z_0_xxyyz_yzzz[k] = -g_z_0_xyyz_yzzz[k] * cd_x[k] + g_z_0_xyyz_xyzzz[k];

                g_z_0_xxyyz_zzzz[k] = -g_z_0_xyyz_zzzz[k] * cd_x[k] + g_z_0_xyyz_xzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 120);

            auto g_z_0_xxyzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 121);

            auto g_z_0_xxyzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 122);

            auto g_z_0_xxyzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 123);

            auto g_z_0_xxyzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 124);

            auto g_z_0_xxyzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 125);

            auto g_z_0_xxyzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 126);

            auto g_z_0_xxyzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 127);

            auto g_z_0_xxyzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 128);

            auto g_z_0_xxyzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 129);

            auto g_z_0_xxyzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 130);

            auto g_z_0_xxyzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 131);

            auto g_z_0_xxyzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 132);

            auto g_z_0_xxyzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 133);

            auto g_z_0_xxyzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 134);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzz_xxxx, g_z_0_xxyzz_xxxy, g_z_0_xxyzz_xxxz, g_z_0_xxyzz_xxyy, g_z_0_xxyzz_xxyz, g_z_0_xxyzz_xxzz, g_z_0_xxyzz_xyyy, g_z_0_xxyzz_xyyz, g_z_0_xxyzz_xyzz, g_z_0_xxyzz_xzzz, g_z_0_xxyzz_yyyy, g_z_0_xxyzz_yyyz, g_z_0_xxyzz_yyzz, g_z_0_xxyzz_yzzz, g_z_0_xxyzz_zzzz, g_z_0_xyzz_xxxx, g_z_0_xyzz_xxxxx, g_z_0_xyzz_xxxxy, g_z_0_xyzz_xxxxz, g_z_0_xyzz_xxxy, g_z_0_xyzz_xxxyy, g_z_0_xyzz_xxxyz, g_z_0_xyzz_xxxz, g_z_0_xyzz_xxxzz, g_z_0_xyzz_xxyy, g_z_0_xyzz_xxyyy, g_z_0_xyzz_xxyyz, g_z_0_xyzz_xxyz, g_z_0_xyzz_xxyzz, g_z_0_xyzz_xxzz, g_z_0_xyzz_xxzzz, g_z_0_xyzz_xyyy, g_z_0_xyzz_xyyyy, g_z_0_xyzz_xyyyz, g_z_0_xyzz_xyyz, g_z_0_xyzz_xyyzz, g_z_0_xyzz_xyzz, g_z_0_xyzz_xyzzz, g_z_0_xyzz_xzzz, g_z_0_xyzz_xzzzz, g_z_0_xyzz_yyyy, g_z_0_xyzz_yyyz, g_z_0_xyzz_yyzz, g_z_0_xyzz_yzzz, g_z_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_xxxx[k] = -g_z_0_xyzz_xxxx[k] * cd_x[k] + g_z_0_xyzz_xxxxx[k];

                g_z_0_xxyzz_xxxy[k] = -g_z_0_xyzz_xxxy[k] * cd_x[k] + g_z_0_xyzz_xxxxy[k];

                g_z_0_xxyzz_xxxz[k] = -g_z_0_xyzz_xxxz[k] * cd_x[k] + g_z_0_xyzz_xxxxz[k];

                g_z_0_xxyzz_xxyy[k] = -g_z_0_xyzz_xxyy[k] * cd_x[k] + g_z_0_xyzz_xxxyy[k];

                g_z_0_xxyzz_xxyz[k] = -g_z_0_xyzz_xxyz[k] * cd_x[k] + g_z_0_xyzz_xxxyz[k];

                g_z_0_xxyzz_xxzz[k] = -g_z_0_xyzz_xxzz[k] * cd_x[k] + g_z_0_xyzz_xxxzz[k];

                g_z_0_xxyzz_xyyy[k] = -g_z_0_xyzz_xyyy[k] * cd_x[k] + g_z_0_xyzz_xxyyy[k];

                g_z_0_xxyzz_xyyz[k] = -g_z_0_xyzz_xyyz[k] * cd_x[k] + g_z_0_xyzz_xxyyz[k];

                g_z_0_xxyzz_xyzz[k] = -g_z_0_xyzz_xyzz[k] * cd_x[k] + g_z_0_xyzz_xxyzz[k];

                g_z_0_xxyzz_xzzz[k] = -g_z_0_xyzz_xzzz[k] * cd_x[k] + g_z_0_xyzz_xxzzz[k];

                g_z_0_xxyzz_yyyy[k] = -g_z_0_xyzz_yyyy[k] * cd_x[k] + g_z_0_xyzz_xyyyy[k];

                g_z_0_xxyzz_yyyz[k] = -g_z_0_xyzz_yyyz[k] * cd_x[k] + g_z_0_xyzz_xyyyz[k];

                g_z_0_xxyzz_yyzz[k] = -g_z_0_xyzz_yyzz[k] * cd_x[k] + g_z_0_xyzz_xyyzz[k];

                g_z_0_xxyzz_yzzz[k] = -g_z_0_xyzz_yzzz[k] * cd_x[k] + g_z_0_xyzz_xyzzz[k];

                g_z_0_xxyzz_zzzz[k] = -g_z_0_xyzz_zzzz[k] * cd_x[k] + g_z_0_xyzz_xzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 135);

            auto g_z_0_xxzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 136);

            auto g_z_0_xxzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 137);

            auto g_z_0_xxzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 138);

            auto g_z_0_xxzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 139);

            auto g_z_0_xxzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 140);

            auto g_z_0_xxzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 141);

            auto g_z_0_xxzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 142);

            auto g_z_0_xxzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 143);

            auto g_z_0_xxzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 144);

            auto g_z_0_xxzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 145);

            auto g_z_0_xxzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 146);

            auto g_z_0_xxzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 147);

            auto g_z_0_xxzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 148);

            auto g_z_0_xxzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzz_xxxx, g_z_0_xxzzz_xxxy, g_z_0_xxzzz_xxxz, g_z_0_xxzzz_xxyy, g_z_0_xxzzz_xxyz, g_z_0_xxzzz_xxzz, g_z_0_xxzzz_xyyy, g_z_0_xxzzz_xyyz, g_z_0_xxzzz_xyzz, g_z_0_xxzzz_xzzz, g_z_0_xxzzz_yyyy, g_z_0_xxzzz_yyyz, g_z_0_xxzzz_yyzz, g_z_0_xxzzz_yzzz, g_z_0_xxzzz_zzzz, g_z_0_xzzz_xxxx, g_z_0_xzzz_xxxxx, g_z_0_xzzz_xxxxy, g_z_0_xzzz_xxxxz, g_z_0_xzzz_xxxy, g_z_0_xzzz_xxxyy, g_z_0_xzzz_xxxyz, g_z_0_xzzz_xxxz, g_z_0_xzzz_xxxzz, g_z_0_xzzz_xxyy, g_z_0_xzzz_xxyyy, g_z_0_xzzz_xxyyz, g_z_0_xzzz_xxyz, g_z_0_xzzz_xxyzz, g_z_0_xzzz_xxzz, g_z_0_xzzz_xxzzz, g_z_0_xzzz_xyyy, g_z_0_xzzz_xyyyy, g_z_0_xzzz_xyyyz, g_z_0_xzzz_xyyz, g_z_0_xzzz_xyyzz, g_z_0_xzzz_xyzz, g_z_0_xzzz_xyzzz, g_z_0_xzzz_xzzz, g_z_0_xzzz_xzzzz, g_z_0_xzzz_yyyy, g_z_0_xzzz_yyyz, g_z_0_xzzz_yyzz, g_z_0_xzzz_yzzz, g_z_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_xxxx[k] = -g_z_0_xzzz_xxxx[k] * cd_x[k] + g_z_0_xzzz_xxxxx[k];

                g_z_0_xxzzz_xxxy[k] = -g_z_0_xzzz_xxxy[k] * cd_x[k] + g_z_0_xzzz_xxxxy[k];

                g_z_0_xxzzz_xxxz[k] = -g_z_0_xzzz_xxxz[k] * cd_x[k] + g_z_0_xzzz_xxxxz[k];

                g_z_0_xxzzz_xxyy[k] = -g_z_0_xzzz_xxyy[k] * cd_x[k] + g_z_0_xzzz_xxxyy[k];

                g_z_0_xxzzz_xxyz[k] = -g_z_0_xzzz_xxyz[k] * cd_x[k] + g_z_0_xzzz_xxxyz[k];

                g_z_0_xxzzz_xxzz[k] = -g_z_0_xzzz_xxzz[k] * cd_x[k] + g_z_0_xzzz_xxxzz[k];

                g_z_0_xxzzz_xyyy[k] = -g_z_0_xzzz_xyyy[k] * cd_x[k] + g_z_0_xzzz_xxyyy[k];

                g_z_0_xxzzz_xyyz[k] = -g_z_0_xzzz_xyyz[k] * cd_x[k] + g_z_0_xzzz_xxyyz[k];

                g_z_0_xxzzz_xyzz[k] = -g_z_0_xzzz_xyzz[k] * cd_x[k] + g_z_0_xzzz_xxyzz[k];

                g_z_0_xxzzz_xzzz[k] = -g_z_0_xzzz_xzzz[k] * cd_x[k] + g_z_0_xzzz_xxzzz[k];

                g_z_0_xxzzz_yyyy[k] = -g_z_0_xzzz_yyyy[k] * cd_x[k] + g_z_0_xzzz_xyyyy[k];

                g_z_0_xxzzz_yyyz[k] = -g_z_0_xzzz_yyyz[k] * cd_x[k] + g_z_0_xzzz_xyyyz[k];

                g_z_0_xxzzz_yyzz[k] = -g_z_0_xzzz_yyzz[k] * cd_x[k] + g_z_0_xzzz_xyyzz[k];

                g_z_0_xxzzz_yzzz[k] = -g_z_0_xzzz_yzzz[k] * cd_x[k] + g_z_0_xzzz_xyzzz[k];

                g_z_0_xxzzz_zzzz[k] = -g_z_0_xzzz_zzzz[k] * cd_x[k] + g_z_0_xzzz_xzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyy_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 150);

            auto g_z_0_xyyyy_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 151);

            auto g_z_0_xyyyy_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 152);

            auto g_z_0_xyyyy_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 153);

            auto g_z_0_xyyyy_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 154);

            auto g_z_0_xyyyy_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 155);

            auto g_z_0_xyyyy_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 156);

            auto g_z_0_xyyyy_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 157);

            auto g_z_0_xyyyy_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 158);

            auto g_z_0_xyyyy_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 159);

            auto g_z_0_xyyyy_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 160);

            auto g_z_0_xyyyy_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 161);

            auto g_z_0_xyyyy_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 162);

            auto g_z_0_xyyyy_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 163);

            auto g_z_0_xyyyy_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 164);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyy_xxxx, g_z_0_xyyyy_xxxy, g_z_0_xyyyy_xxxz, g_z_0_xyyyy_xxyy, g_z_0_xyyyy_xxyz, g_z_0_xyyyy_xxzz, g_z_0_xyyyy_xyyy, g_z_0_xyyyy_xyyz, g_z_0_xyyyy_xyzz, g_z_0_xyyyy_xzzz, g_z_0_xyyyy_yyyy, g_z_0_xyyyy_yyyz, g_z_0_xyyyy_yyzz, g_z_0_xyyyy_yzzz, g_z_0_xyyyy_zzzz, g_z_0_yyyy_xxxx, g_z_0_yyyy_xxxxx, g_z_0_yyyy_xxxxy, g_z_0_yyyy_xxxxz, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxxyy, g_z_0_yyyy_xxxyz, g_z_0_yyyy_xxxz, g_z_0_yyyy_xxxzz, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyyy, g_z_0_yyyy_xxyyz, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxyzz, g_z_0_yyyy_xxzz, g_z_0_yyyy_xxzzz, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyyy, g_z_0_yyyy_xyyyz, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyyzz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xyzzz, g_z_0_yyyy_xzzz, g_z_0_yyyy_xzzzz, g_z_0_yyyy_yyyy, g_z_0_yyyy_yyyz, g_z_0_yyyy_yyzz, g_z_0_yyyy_yzzz, g_z_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_xxxx[k] = -g_z_0_yyyy_xxxx[k] * cd_x[k] + g_z_0_yyyy_xxxxx[k];

                g_z_0_xyyyy_xxxy[k] = -g_z_0_yyyy_xxxy[k] * cd_x[k] + g_z_0_yyyy_xxxxy[k];

                g_z_0_xyyyy_xxxz[k] = -g_z_0_yyyy_xxxz[k] * cd_x[k] + g_z_0_yyyy_xxxxz[k];

                g_z_0_xyyyy_xxyy[k] = -g_z_0_yyyy_xxyy[k] * cd_x[k] + g_z_0_yyyy_xxxyy[k];

                g_z_0_xyyyy_xxyz[k] = -g_z_0_yyyy_xxyz[k] * cd_x[k] + g_z_0_yyyy_xxxyz[k];

                g_z_0_xyyyy_xxzz[k] = -g_z_0_yyyy_xxzz[k] * cd_x[k] + g_z_0_yyyy_xxxzz[k];

                g_z_0_xyyyy_xyyy[k] = -g_z_0_yyyy_xyyy[k] * cd_x[k] + g_z_0_yyyy_xxyyy[k];

                g_z_0_xyyyy_xyyz[k] = -g_z_0_yyyy_xyyz[k] * cd_x[k] + g_z_0_yyyy_xxyyz[k];

                g_z_0_xyyyy_xyzz[k] = -g_z_0_yyyy_xyzz[k] * cd_x[k] + g_z_0_yyyy_xxyzz[k];

                g_z_0_xyyyy_xzzz[k] = -g_z_0_yyyy_xzzz[k] * cd_x[k] + g_z_0_yyyy_xxzzz[k];

                g_z_0_xyyyy_yyyy[k] = -g_z_0_yyyy_yyyy[k] * cd_x[k] + g_z_0_yyyy_xyyyy[k];

                g_z_0_xyyyy_yyyz[k] = -g_z_0_yyyy_yyyz[k] * cd_x[k] + g_z_0_yyyy_xyyyz[k];

                g_z_0_xyyyy_yyzz[k] = -g_z_0_yyyy_yyzz[k] * cd_x[k] + g_z_0_yyyy_xyyzz[k];

                g_z_0_xyyyy_yzzz[k] = -g_z_0_yyyy_yzzz[k] * cd_x[k] + g_z_0_yyyy_xyzzz[k];

                g_z_0_xyyyy_zzzz[k] = -g_z_0_yyyy_zzzz[k] * cd_x[k] + g_z_0_yyyy_xzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 165);

            auto g_z_0_xyyyz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 166);

            auto g_z_0_xyyyz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 167);

            auto g_z_0_xyyyz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 168);

            auto g_z_0_xyyyz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 169);

            auto g_z_0_xyyyz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 170);

            auto g_z_0_xyyyz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 171);

            auto g_z_0_xyyyz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 172);

            auto g_z_0_xyyyz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 173);

            auto g_z_0_xyyyz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 174);

            auto g_z_0_xyyyz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 175);

            auto g_z_0_xyyyz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 176);

            auto g_z_0_xyyyz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 177);

            auto g_z_0_xyyyz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 178);

            auto g_z_0_xyyyz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyz_xxxx, g_z_0_xyyyz_xxxy, g_z_0_xyyyz_xxxz, g_z_0_xyyyz_xxyy, g_z_0_xyyyz_xxyz, g_z_0_xyyyz_xxzz, g_z_0_xyyyz_xyyy, g_z_0_xyyyz_xyyz, g_z_0_xyyyz_xyzz, g_z_0_xyyyz_xzzz, g_z_0_xyyyz_yyyy, g_z_0_xyyyz_yyyz, g_z_0_xyyyz_yyzz, g_z_0_xyyyz_yzzz, g_z_0_xyyyz_zzzz, g_z_0_yyyz_xxxx, g_z_0_yyyz_xxxxx, g_z_0_yyyz_xxxxy, g_z_0_yyyz_xxxxz, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxxyy, g_z_0_yyyz_xxxyz, g_z_0_yyyz_xxxz, g_z_0_yyyz_xxxzz, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyyy, g_z_0_yyyz_xxyyz, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxyzz, g_z_0_yyyz_xxzz, g_z_0_yyyz_xxzzz, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyyy, g_z_0_yyyz_xyyyz, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyyzz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xyzzz, g_z_0_yyyz_xzzz, g_z_0_yyyz_xzzzz, g_z_0_yyyz_yyyy, g_z_0_yyyz_yyyz, g_z_0_yyyz_yyzz, g_z_0_yyyz_yzzz, g_z_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_xxxx[k] = -g_z_0_yyyz_xxxx[k] * cd_x[k] + g_z_0_yyyz_xxxxx[k];

                g_z_0_xyyyz_xxxy[k] = -g_z_0_yyyz_xxxy[k] * cd_x[k] + g_z_0_yyyz_xxxxy[k];

                g_z_0_xyyyz_xxxz[k] = -g_z_0_yyyz_xxxz[k] * cd_x[k] + g_z_0_yyyz_xxxxz[k];

                g_z_0_xyyyz_xxyy[k] = -g_z_0_yyyz_xxyy[k] * cd_x[k] + g_z_0_yyyz_xxxyy[k];

                g_z_0_xyyyz_xxyz[k] = -g_z_0_yyyz_xxyz[k] * cd_x[k] + g_z_0_yyyz_xxxyz[k];

                g_z_0_xyyyz_xxzz[k] = -g_z_0_yyyz_xxzz[k] * cd_x[k] + g_z_0_yyyz_xxxzz[k];

                g_z_0_xyyyz_xyyy[k] = -g_z_0_yyyz_xyyy[k] * cd_x[k] + g_z_0_yyyz_xxyyy[k];

                g_z_0_xyyyz_xyyz[k] = -g_z_0_yyyz_xyyz[k] * cd_x[k] + g_z_0_yyyz_xxyyz[k];

                g_z_0_xyyyz_xyzz[k] = -g_z_0_yyyz_xyzz[k] * cd_x[k] + g_z_0_yyyz_xxyzz[k];

                g_z_0_xyyyz_xzzz[k] = -g_z_0_yyyz_xzzz[k] * cd_x[k] + g_z_0_yyyz_xxzzz[k];

                g_z_0_xyyyz_yyyy[k] = -g_z_0_yyyz_yyyy[k] * cd_x[k] + g_z_0_yyyz_xyyyy[k];

                g_z_0_xyyyz_yyyz[k] = -g_z_0_yyyz_yyyz[k] * cd_x[k] + g_z_0_yyyz_xyyyz[k];

                g_z_0_xyyyz_yyzz[k] = -g_z_0_yyyz_yyzz[k] * cd_x[k] + g_z_0_yyyz_xyyzz[k];

                g_z_0_xyyyz_yzzz[k] = -g_z_0_yyyz_yzzz[k] * cd_x[k] + g_z_0_yyyz_xyzzz[k];

                g_z_0_xyyyz_zzzz[k] = -g_z_0_yyyz_zzzz[k] * cd_x[k] + g_z_0_yyyz_xzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 180);

            auto g_z_0_xyyzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 181);

            auto g_z_0_xyyzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 182);

            auto g_z_0_xyyzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 183);

            auto g_z_0_xyyzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 184);

            auto g_z_0_xyyzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 185);

            auto g_z_0_xyyzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 186);

            auto g_z_0_xyyzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 187);

            auto g_z_0_xyyzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 188);

            auto g_z_0_xyyzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 189);

            auto g_z_0_xyyzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 190);

            auto g_z_0_xyyzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 191);

            auto g_z_0_xyyzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 192);

            auto g_z_0_xyyzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 193);

            auto g_z_0_xyyzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 194);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzz_xxxx, g_z_0_xyyzz_xxxy, g_z_0_xyyzz_xxxz, g_z_0_xyyzz_xxyy, g_z_0_xyyzz_xxyz, g_z_0_xyyzz_xxzz, g_z_0_xyyzz_xyyy, g_z_0_xyyzz_xyyz, g_z_0_xyyzz_xyzz, g_z_0_xyyzz_xzzz, g_z_0_xyyzz_yyyy, g_z_0_xyyzz_yyyz, g_z_0_xyyzz_yyzz, g_z_0_xyyzz_yzzz, g_z_0_xyyzz_zzzz, g_z_0_yyzz_xxxx, g_z_0_yyzz_xxxxx, g_z_0_yyzz_xxxxy, g_z_0_yyzz_xxxxz, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxxyy, g_z_0_yyzz_xxxyz, g_z_0_yyzz_xxxz, g_z_0_yyzz_xxxzz, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyyy, g_z_0_yyzz_xxyyz, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxyzz, g_z_0_yyzz_xxzz, g_z_0_yyzz_xxzzz, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyyy, g_z_0_yyzz_xyyyz, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyyzz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xyzzz, g_z_0_yyzz_xzzz, g_z_0_yyzz_xzzzz, g_z_0_yyzz_yyyy, g_z_0_yyzz_yyyz, g_z_0_yyzz_yyzz, g_z_0_yyzz_yzzz, g_z_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_xxxx[k] = -g_z_0_yyzz_xxxx[k] * cd_x[k] + g_z_0_yyzz_xxxxx[k];

                g_z_0_xyyzz_xxxy[k] = -g_z_0_yyzz_xxxy[k] * cd_x[k] + g_z_0_yyzz_xxxxy[k];

                g_z_0_xyyzz_xxxz[k] = -g_z_0_yyzz_xxxz[k] * cd_x[k] + g_z_0_yyzz_xxxxz[k];

                g_z_0_xyyzz_xxyy[k] = -g_z_0_yyzz_xxyy[k] * cd_x[k] + g_z_0_yyzz_xxxyy[k];

                g_z_0_xyyzz_xxyz[k] = -g_z_0_yyzz_xxyz[k] * cd_x[k] + g_z_0_yyzz_xxxyz[k];

                g_z_0_xyyzz_xxzz[k] = -g_z_0_yyzz_xxzz[k] * cd_x[k] + g_z_0_yyzz_xxxzz[k];

                g_z_0_xyyzz_xyyy[k] = -g_z_0_yyzz_xyyy[k] * cd_x[k] + g_z_0_yyzz_xxyyy[k];

                g_z_0_xyyzz_xyyz[k] = -g_z_0_yyzz_xyyz[k] * cd_x[k] + g_z_0_yyzz_xxyyz[k];

                g_z_0_xyyzz_xyzz[k] = -g_z_0_yyzz_xyzz[k] * cd_x[k] + g_z_0_yyzz_xxyzz[k];

                g_z_0_xyyzz_xzzz[k] = -g_z_0_yyzz_xzzz[k] * cd_x[k] + g_z_0_yyzz_xxzzz[k];

                g_z_0_xyyzz_yyyy[k] = -g_z_0_yyzz_yyyy[k] * cd_x[k] + g_z_0_yyzz_xyyyy[k];

                g_z_0_xyyzz_yyyz[k] = -g_z_0_yyzz_yyyz[k] * cd_x[k] + g_z_0_yyzz_xyyyz[k];

                g_z_0_xyyzz_yyzz[k] = -g_z_0_yyzz_yyzz[k] * cd_x[k] + g_z_0_yyzz_xyyzz[k];

                g_z_0_xyyzz_yzzz[k] = -g_z_0_yyzz_yzzz[k] * cd_x[k] + g_z_0_yyzz_xyzzz[k];

                g_z_0_xyyzz_zzzz[k] = -g_z_0_yyzz_zzzz[k] * cd_x[k] + g_z_0_yyzz_xzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 195);

            auto g_z_0_xyzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 196);

            auto g_z_0_xyzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 197);

            auto g_z_0_xyzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 198);

            auto g_z_0_xyzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 199);

            auto g_z_0_xyzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 200);

            auto g_z_0_xyzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 201);

            auto g_z_0_xyzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 202);

            auto g_z_0_xyzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 203);

            auto g_z_0_xyzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 204);

            auto g_z_0_xyzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 205);

            auto g_z_0_xyzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 206);

            auto g_z_0_xyzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 207);

            auto g_z_0_xyzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 208);

            auto g_z_0_xyzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzz_xxxx, g_z_0_xyzzz_xxxy, g_z_0_xyzzz_xxxz, g_z_0_xyzzz_xxyy, g_z_0_xyzzz_xxyz, g_z_0_xyzzz_xxzz, g_z_0_xyzzz_xyyy, g_z_0_xyzzz_xyyz, g_z_0_xyzzz_xyzz, g_z_0_xyzzz_xzzz, g_z_0_xyzzz_yyyy, g_z_0_xyzzz_yyyz, g_z_0_xyzzz_yyzz, g_z_0_xyzzz_yzzz, g_z_0_xyzzz_zzzz, g_z_0_yzzz_xxxx, g_z_0_yzzz_xxxxx, g_z_0_yzzz_xxxxy, g_z_0_yzzz_xxxxz, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxxyy, g_z_0_yzzz_xxxyz, g_z_0_yzzz_xxxz, g_z_0_yzzz_xxxzz, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyyy, g_z_0_yzzz_xxyyz, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxyzz, g_z_0_yzzz_xxzz, g_z_0_yzzz_xxzzz, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyyy, g_z_0_yzzz_xyyyz, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyyzz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xyzzz, g_z_0_yzzz_xzzz, g_z_0_yzzz_xzzzz, g_z_0_yzzz_yyyy, g_z_0_yzzz_yyyz, g_z_0_yzzz_yyzz, g_z_0_yzzz_yzzz, g_z_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_xxxx[k] = -g_z_0_yzzz_xxxx[k] * cd_x[k] + g_z_0_yzzz_xxxxx[k];

                g_z_0_xyzzz_xxxy[k] = -g_z_0_yzzz_xxxy[k] * cd_x[k] + g_z_0_yzzz_xxxxy[k];

                g_z_0_xyzzz_xxxz[k] = -g_z_0_yzzz_xxxz[k] * cd_x[k] + g_z_0_yzzz_xxxxz[k];

                g_z_0_xyzzz_xxyy[k] = -g_z_0_yzzz_xxyy[k] * cd_x[k] + g_z_0_yzzz_xxxyy[k];

                g_z_0_xyzzz_xxyz[k] = -g_z_0_yzzz_xxyz[k] * cd_x[k] + g_z_0_yzzz_xxxyz[k];

                g_z_0_xyzzz_xxzz[k] = -g_z_0_yzzz_xxzz[k] * cd_x[k] + g_z_0_yzzz_xxxzz[k];

                g_z_0_xyzzz_xyyy[k] = -g_z_0_yzzz_xyyy[k] * cd_x[k] + g_z_0_yzzz_xxyyy[k];

                g_z_0_xyzzz_xyyz[k] = -g_z_0_yzzz_xyyz[k] * cd_x[k] + g_z_0_yzzz_xxyyz[k];

                g_z_0_xyzzz_xyzz[k] = -g_z_0_yzzz_xyzz[k] * cd_x[k] + g_z_0_yzzz_xxyzz[k];

                g_z_0_xyzzz_xzzz[k] = -g_z_0_yzzz_xzzz[k] * cd_x[k] + g_z_0_yzzz_xxzzz[k];

                g_z_0_xyzzz_yyyy[k] = -g_z_0_yzzz_yyyy[k] * cd_x[k] + g_z_0_yzzz_xyyyy[k];

                g_z_0_xyzzz_yyyz[k] = -g_z_0_yzzz_yyyz[k] * cd_x[k] + g_z_0_yzzz_xyyyz[k];

                g_z_0_xyzzz_yyzz[k] = -g_z_0_yzzz_yyzz[k] * cd_x[k] + g_z_0_yzzz_xyyzz[k];

                g_z_0_xyzzz_yzzz[k] = -g_z_0_yzzz_yzzz[k] * cd_x[k] + g_z_0_yzzz_xyzzz[k];

                g_z_0_xyzzz_zzzz[k] = -g_z_0_yzzz_zzzz[k] * cd_x[k] + g_z_0_yzzz_xzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 210);

            auto g_z_0_xzzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 211);

            auto g_z_0_xzzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 212);

            auto g_z_0_xzzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 213);

            auto g_z_0_xzzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 214);

            auto g_z_0_xzzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 215);

            auto g_z_0_xzzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 216);

            auto g_z_0_xzzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 217);

            auto g_z_0_xzzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 218);

            auto g_z_0_xzzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 219);

            auto g_z_0_xzzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 220);

            auto g_z_0_xzzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 221);

            auto g_z_0_xzzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 222);

            auto g_z_0_xzzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 223);

            auto g_z_0_xzzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 224);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzz_xxxx, g_z_0_xzzzz_xxxy, g_z_0_xzzzz_xxxz, g_z_0_xzzzz_xxyy, g_z_0_xzzzz_xxyz, g_z_0_xzzzz_xxzz, g_z_0_xzzzz_xyyy, g_z_0_xzzzz_xyyz, g_z_0_xzzzz_xyzz, g_z_0_xzzzz_xzzz, g_z_0_xzzzz_yyyy, g_z_0_xzzzz_yyyz, g_z_0_xzzzz_yyzz, g_z_0_xzzzz_yzzz, g_z_0_xzzzz_zzzz, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxxx, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_xxxx[k] = -g_z_0_zzzz_xxxx[k] * cd_x[k] + g_z_0_zzzz_xxxxx[k];

                g_z_0_xzzzz_xxxy[k] = -g_z_0_zzzz_xxxy[k] * cd_x[k] + g_z_0_zzzz_xxxxy[k];

                g_z_0_xzzzz_xxxz[k] = -g_z_0_zzzz_xxxz[k] * cd_x[k] + g_z_0_zzzz_xxxxz[k];

                g_z_0_xzzzz_xxyy[k] = -g_z_0_zzzz_xxyy[k] * cd_x[k] + g_z_0_zzzz_xxxyy[k];

                g_z_0_xzzzz_xxyz[k] = -g_z_0_zzzz_xxyz[k] * cd_x[k] + g_z_0_zzzz_xxxyz[k];

                g_z_0_xzzzz_xxzz[k] = -g_z_0_zzzz_xxzz[k] * cd_x[k] + g_z_0_zzzz_xxxzz[k];

                g_z_0_xzzzz_xyyy[k] = -g_z_0_zzzz_xyyy[k] * cd_x[k] + g_z_0_zzzz_xxyyy[k];

                g_z_0_xzzzz_xyyz[k] = -g_z_0_zzzz_xyyz[k] * cd_x[k] + g_z_0_zzzz_xxyyz[k];

                g_z_0_xzzzz_xyzz[k] = -g_z_0_zzzz_xyzz[k] * cd_x[k] + g_z_0_zzzz_xxyzz[k];

                g_z_0_xzzzz_xzzz[k] = -g_z_0_zzzz_xzzz[k] * cd_x[k] + g_z_0_zzzz_xxzzz[k];

                g_z_0_xzzzz_yyyy[k] = -g_z_0_zzzz_yyyy[k] * cd_x[k] + g_z_0_zzzz_xyyyy[k];

                g_z_0_xzzzz_yyyz[k] = -g_z_0_zzzz_yyyz[k] * cd_x[k] + g_z_0_zzzz_xyyyz[k];

                g_z_0_xzzzz_yyzz[k] = -g_z_0_zzzz_yyzz[k] * cd_x[k] + g_z_0_zzzz_xyyzz[k];

                g_z_0_xzzzz_yzzz[k] = -g_z_0_zzzz_yzzz[k] * cd_x[k] + g_z_0_zzzz_xyzzz[k];

                g_z_0_xzzzz_zzzz[k] = -g_z_0_zzzz_zzzz[k] * cd_x[k] + g_z_0_zzzz_xzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyy_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 225);

            auto g_z_0_yyyyy_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 226);

            auto g_z_0_yyyyy_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 227);

            auto g_z_0_yyyyy_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 228);

            auto g_z_0_yyyyy_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 229);

            auto g_z_0_yyyyy_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 230);

            auto g_z_0_yyyyy_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 231);

            auto g_z_0_yyyyy_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 232);

            auto g_z_0_yyyyy_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 233);

            auto g_z_0_yyyyy_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 234);

            auto g_z_0_yyyyy_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 235);

            auto g_z_0_yyyyy_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 236);

            auto g_z_0_yyyyy_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 237);

            auto g_z_0_yyyyy_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 238);

            auto g_z_0_yyyyy_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 239);

            #pragma omp simd aligned(cd_y, g_z_0_yyyy_xxxx, g_z_0_yyyy_xxxxy, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxxyy, g_z_0_yyyy_xxxyz, g_z_0_yyyy_xxxz, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyyy, g_z_0_yyyy_xxyyz, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxyzz, g_z_0_yyyy_xxzz, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyyy, g_z_0_yyyy_xyyyz, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyyzz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xyzzz, g_z_0_yyyy_xzzz, g_z_0_yyyy_yyyy, g_z_0_yyyy_yyyyy, g_z_0_yyyy_yyyyz, g_z_0_yyyy_yyyz, g_z_0_yyyy_yyyzz, g_z_0_yyyy_yyzz, g_z_0_yyyy_yyzzz, g_z_0_yyyy_yzzz, g_z_0_yyyy_yzzzz, g_z_0_yyyy_zzzz, g_z_0_yyyyy_xxxx, g_z_0_yyyyy_xxxy, g_z_0_yyyyy_xxxz, g_z_0_yyyyy_xxyy, g_z_0_yyyyy_xxyz, g_z_0_yyyyy_xxzz, g_z_0_yyyyy_xyyy, g_z_0_yyyyy_xyyz, g_z_0_yyyyy_xyzz, g_z_0_yyyyy_xzzz, g_z_0_yyyyy_yyyy, g_z_0_yyyyy_yyyz, g_z_0_yyyyy_yyzz, g_z_0_yyyyy_yzzz, g_z_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_xxxx[k] = -g_z_0_yyyy_xxxx[k] * cd_y[k] + g_z_0_yyyy_xxxxy[k];

                g_z_0_yyyyy_xxxy[k] = -g_z_0_yyyy_xxxy[k] * cd_y[k] + g_z_0_yyyy_xxxyy[k];

                g_z_0_yyyyy_xxxz[k] = -g_z_0_yyyy_xxxz[k] * cd_y[k] + g_z_0_yyyy_xxxyz[k];

                g_z_0_yyyyy_xxyy[k] = -g_z_0_yyyy_xxyy[k] * cd_y[k] + g_z_0_yyyy_xxyyy[k];

                g_z_0_yyyyy_xxyz[k] = -g_z_0_yyyy_xxyz[k] * cd_y[k] + g_z_0_yyyy_xxyyz[k];

                g_z_0_yyyyy_xxzz[k] = -g_z_0_yyyy_xxzz[k] * cd_y[k] + g_z_0_yyyy_xxyzz[k];

                g_z_0_yyyyy_xyyy[k] = -g_z_0_yyyy_xyyy[k] * cd_y[k] + g_z_0_yyyy_xyyyy[k];

                g_z_0_yyyyy_xyyz[k] = -g_z_0_yyyy_xyyz[k] * cd_y[k] + g_z_0_yyyy_xyyyz[k];

                g_z_0_yyyyy_xyzz[k] = -g_z_0_yyyy_xyzz[k] * cd_y[k] + g_z_0_yyyy_xyyzz[k];

                g_z_0_yyyyy_xzzz[k] = -g_z_0_yyyy_xzzz[k] * cd_y[k] + g_z_0_yyyy_xyzzz[k];

                g_z_0_yyyyy_yyyy[k] = -g_z_0_yyyy_yyyy[k] * cd_y[k] + g_z_0_yyyy_yyyyy[k];

                g_z_0_yyyyy_yyyz[k] = -g_z_0_yyyy_yyyz[k] * cd_y[k] + g_z_0_yyyy_yyyyz[k];

                g_z_0_yyyyy_yyzz[k] = -g_z_0_yyyy_yyzz[k] * cd_y[k] + g_z_0_yyyy_yyyzz[k];

                g_z_0_yyyyy_yzzz[k] = -g_z_0_yyyy_yzzz[k] * cd_y[k] + g_z_0_yyyy_yyzzz[k];

                g_z_0_yyyyy_zzzz[k] = -g_z_0_yyyy_zzzz[k] * cd_y[k] + g_z_0_yyyy_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 240);

            auto g_z_0_yyyyz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 241);

            auto g_z_0_yyyyz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 242);

            auto g_z_0_yyyyz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 243);

            auto g_z_0_yyyyz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 244);

            auto g_z_0_yyyyz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 245);

            auto g_z_0_yyyyz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 246);

            auto g_z_0_yyyyz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 247);

            auto g_z_0_yyyyz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 248);

            auto g_z_0_yyyyz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 249);

            auto g_z_0_yyyyz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 250);

            auto g_z_0_yyyyz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 251);

            auto g_z_0_yyyyz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 252);

            auto g_z_0_yyyyz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 253);

            auto g_z_0_yyyyz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 254);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyz_xxxx, g_z_0_yyyyz_xxxy, g_z_0_yyyyz_xxxz, g_z_0_yyyyz_xxyy, g_z_0_yyyyz_xxyz, g_z_0_yyyyz_xxzz, g_z_0_yyyyz_xyyy, g_z_0_yyyyz_xyyz, g_z_0_yyyyz_xyzz, g_z_0_yyyyz_xzzz, g_z_0_yyyyz_yyyy, g_z_0_yyyyz_yyyz, g_z_0_yyyyz_yyzz, g_z_0_yyyyz_yzzz, g_z_0_yyyyz_zzzz, g_z_0_yyyz_xxxx, g_z_0_yyyz_xxxxy, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxxyy, g_z_0_yyyz_xxxyz, g_z_0_yyyz_xxxz, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyyy, g_z_0_yyyz_xxyyz, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxyzz, g_z_0_yyyz_xxzz, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyyy, g_z_0_yyyz_xyyyz, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyyzz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xyzzz, g_z_0_yyyz_xzzz, g_z_0_yyyz_yyyy, g_z_0_yyyz_yyyyy, g_z_0_yyyz_yyyyz, g_z_0_yyyz_yyyz, g_z_0_yyyz_yyyzz, g_z_0_yyyz_yyzz, g_z_0_yyyz_yyzzz, g_z_0_yyyz_yzzz, g_z_0_yyyz_yzzzz, g_z_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_xxxx[k] = -g_z_0_yyyz_xxxx[k] * cd_y[k] + g_z_0_yyyz_xxxxy[k];

                g_z_0_yyyyz_xxxy[k] = -g_z_0_yyyz_xxxy[k] * cd_y[k] + g_z_0_yyyz_xxxyy[k];

                g_z_0_yyyyz_xxxz[k] = -g_z_0_yyyz_xxxz[k] * cd_y[k] + g_z_0_yyyz_xxxyz[k];

                g_z_0_yyyyz_xxyy[k] = -g_z_0_yyyz_xxyy[k] * cd_y[k] + g_z_0_yyyz_xxyyy[k];

                g_z_0_yyyyz_xxyz[k] = -g_z_0_yyyz_xxyz[k] * cd_y[k] + g_z_0_yyyz_xxyyz[k];

                g_z_0_yyyyz_xxzz[k] = -g_z_0_yyyz_xxzz[k] * cd_y[k] + g_z_0_yyyz_xxyzz[k];

                g_z_0_yyyyz_xyyy[k] = -g_z_0_yyyz_xyyy[k] * cd_y[k] + g_z_0_yyyz_xyyyy[k];

                g_z_0_yyyyz_xyyz[k] = -g_z_0_yyyz_xyyz[k] * cd_y[k] + g_z_0_yyyz_xyyyz[k];

                g_z_0_yyyyz_xyzz[k] = -g_z_0_yyyz_xyzz[k] * cd_y[k] + g_z_0_yyyz_xyyzz[k];

                g_z_0_yyyyz_xzzz[k] = -g_z_0_yyyz_xzzz[k] * cd_y[k] + g_z_0_yyyz_xyzzz[k];

                g_z_0_yyyyz_yyyy[k] = -g_z_0_yyyz_yyyy[k] * cd_y[k] + g_z_0_yyyz_yyyyy[k];

                g_z_0_yyyyz_yyyz[k] = -g_z_0_yyyz_yyyz[k] * cd_y[k] + g_z_0_yyyz_yyyyz[k];

                g_z_0_yyyyz_yyzz[k] = -g_z_0_yyyz_yyzz[k] * cd_y[k] + g_z_0_yyyz_yyyzz[k];

                g_z_0_yyyyz_yzzz[k] = -g_z_0_yyyz_yzzz[k] * cd_y[k] + g_z_0_yyyz_yyzzz[k];

                g_z_0_yyyyz_zzzz[k] = -g_z_0_yyyz_zzzz[k] * cd_y[k] + g_z_0_yyyz_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 255);

            auto g_z_0_yyyzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 256);

            auto g_z_0_yyyzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 257);

            auto g_z_0_yyyzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 258);

            auto g_z_0_yyyzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 259);

            auto g_z_0_yyyzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 260);

            auto g_z_0_yyyzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 261);

            auto g_z_0_yyyzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 262);

            auto g_z_0_yyyzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 263);

            auto g_z_0_yyyzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 264);

            auto g_z_0_yyyzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 265);

            auto g_z_0_yyyzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 266);

            auto g_z_0_yyyzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 267);

            auto g_z_0_yyyzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 268);

            auto g_z_0_yyyzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 269);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzz_xxxx, g_z_0_yyyzz_xxxy, g_z_0_yyyzz_xxxz, g_z_0_yyyzz_xxyy, g_z_0_yyyzz_xxyz, g_z_0_yyyzz_xxzz, g_z_0_yyyzz_xyyy, g_z_0_yyyzz_xyyz, g_z_0_yyyzz_xyzz, g_z_0_yyyzz_xzzz, g_z_0_yyyzz_yyyy, g_z_0_yyyzz_yyyz, g_z_0_yyyzz_yyzz, g_z_0_yyyzz_yzzz, g_z_0_yyyzz_zzzz, g_z_0_yyzz_xxxx, g_z_0_yyzz_xxxxy, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxxyy, g_z_0_yyzz_xxxyz, g_z_0_yyzz_xxxz, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyyy, g_z_0_yyzz_xxyyz, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxyzz, g_z_0_yyzz_xxzz, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyyy, g_z_0_yyzz_xyyyz, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyyzz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xyzzz, g_z_0_yyzz_xzzz, g_z_0_yyzz_yyyy, g_z_0_yyzz_yyyyy, g_z_0_yyzz_yyyyz, g_z_0_yyzz_yyyz, g_z_0_yyzz_yyyzz, g_z_0_yyzz_yyzz, g_z_0_yyzz_yyzzz, g_z_0_yyzz_yzzz, g_z_0_yyzz_yzzzz, g_z_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_xxxx[k] = -g_z_0_yyzz_xxxx[k] * cd_y[k] + g_z_0_yyzz_xxxxy[k];

                g_z_0_yyyzz_xxxy[k] = -g_z_0_yyzz_xxxy[k] * cd_y[k] + g_z_0_yyzz_xxxyy[k];

                g_z_0_yyyzz_xxxz[k] = -g_z_0_yyzz_xxxz[k] * cd_y[k] + g_z_0_yyzz_xxxyz[k];

                g_z_0_yyyzz_xxyy[k] = -g_z_0_yyzz_xxyy[k] * cd_y[k] + g_z_0_yyzz_xxyyy[k];

                g_z_0_yyyzz_xxyz[k] = -g_z_0_yyzz_xxyz[k] * cd_y[k] + g_z_0_yyzz_xxyyz[k];

                g_z_0_yyyzz_xxzz[k] = -g_z_0_yyzz_xxzz[k] * cd_y[k] + g_z_0_yyzz_xxyzz[k];

                g_z_0_yyyzz_xyyy[k] = -g_z_0_yyzz_xyyy[k] * cd_y[k] + g_z_0_yyzz_xyyyy[k];

                g_z_0_yyyzz_xyyz[k] = -g_z_0_yyzz_xyyz[k] * cd_y[k] + g_z_0_yyzz_xyyyz[k];

                g_z_0_yyyzz_xyzz[k] = -g_z_0_yyzz_xyzz[k] * cd_y[k] + g_z_0_yyzz_xyyzz[k];

                g_z_0_yyyzz_xzzz[k] = -g_z_0_yyzz_xzzz[k] * cd_y[k] + g_z_0_yyzz_xyzzz[k];

                g_z_0_yyyzz_yyyy[k] = -g_z_0_yyzz_yyyy[k] * cd_y[k] + g_z_0_yyzz_yyyyy[k];

                g_z_0_yyyzz_yyyz[k] = -g_z_0_yyzz_yyyz[k] * cd_y[k] + g_z_0_yyzz_yyyyz[k];

                g_z_0_yyyzz_yyzz[k] = -g_z_0_yyzz_yyzz[k] * cd_y[k] + g_z_0_yyzz_yyyzz[k];

                g_z_0_yyyzz_yzzz[k] = -g_z_0_yyzz_yzzz[k] * cd_y[k] + g_z_0_yyzz_yyzzz[k];

                g_z_0_yyyzz_zzzz[k] = -g_z_0_yyzz_zzzz[k] * cd_y[k] + g_z_0_yyzz_yzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 270);

            auto g_z_0_yyzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 271);

            auto g_z_0_yyzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 272);

            auto g_z_0_yyzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 273);

            auto g_z_0_yyzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 274);

            auto g_z_0_yyzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 275);

            auto g_z_0_yyzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 276);

            auto g_z_0_yyzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 277);

            auto g_z_0_yyzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 278);

            auto g_z_0_yyzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 279);

            auto g_z_0_yyzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 280);

            auto g_z_0_yyzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 281);

            auto g_z_0_yyzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 282);

            auto g_z_0_yyzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 283);

            auto g_z_0_yyzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 284);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzz_xxxx, g_z_0_yyzzz_xxxy, g_z_0_yyzzz_xxxz, g_z_0_yyzzz_xxyy, g_z_0_yyzzz_xxyz, g_z_0_yyzzz_xxzz, g_z_0_yyzzz_xyyy, g_z_0_yyzzz_xyyz, g_z_0_yyzzz_xyzz, g_z_0_yyzzz_xzzz, g_z_0_yyzzz_yyyy, g_z_0_yyzzz_yyyz, g_z_0_yyzzz_yyzz, g_z_0_yyzzz_yzzz, g_z_0_yyzzz_zzzz, g_z_0_yzzz_xxxx, g_z_0_yzzz_xxxxy, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxxyy, g_z_0_yzzz_xxxyz, g_z_0_yzzz_xxxz, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyyy, g_z_0_yzzz_xxyyz, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxyzz, g_z_0_yzzz_xxzz, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyyy, g_z_0_yzzz_xyyyz, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyyzz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xyzzz, g_z_0_yzzz_xzzz, g_z_0_yzzz_yyyy, g_z_0_yzzz_yyyyy, g_z_0_yzzz_yyyyz, g_z_0_yzzz_yyyz, g_z_0_yzzz_yyyzz, g_z_0_yzzz_yyzz, g_z_0_yzzz_yyzzz, g_z_0_yzzz_yzzz, g_z_0_yzzz_yzzzz, g_z_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_xxxx[k] = -g_z_0_yzzz_xxxx[k] * cd_y[k] + g_z_0_yzzz_xxxxy[k];

                g_z_0_yyzzz_xxxy[k] = -g_z_0_yzzz_xxxy[k] * cd_y[k] + g_z_0_yzzz_xxxyy[k];

                g_z_0_yyzzz_xxxz[k] = -g_z_0_yzzz_xxxz[k] * cd_y[k] + g_z_0_yzzz_xxxyz[k];

                g_z_0_yyzzz_xxyy[k] = -g_z_0_yzzz_xxyy[k] * cd_y[k] + g_z_0_yzzz_xxyyy[k];

                g_z_0_yyzzz_xxyz[k] = -g_z_0_yzzz_xxyz[k] * cd_y[k] + g_z_0_yzzz_xxyyz[k];

                g_z_0_yyzzz_xxzz[k] = -g_z_0_yzzz_xxzz[k] * cd_y[k] + g_z_0_yzzz_xxyzz[k];

                g_z_0_yyzzz_xyyy[k] = -g_z_0_yzzz_xyyy[k] * cd_y[k] + g_z_0_yzzz_xyyyy[k];

                g_z_0_yyzzz_xyyz[k] = -g_z_0_yzzz_xyyz[k] * cd_y[k] + g_z_0_yzzz_xyyyz[k];

                g_z_0_yyzzz_xyzz[k] = -g_z_0_yzzz_xyzz[k] * cd_y[k] + g_z_0_yzzz_xyyzz[k];

                g_z_0_yyzzz_xzzz[k] = -g_z_0_yzzz_xzzz[k] * cd_y[k] + g_z_0_yzzz_xyzzz[k];

                g_z_0_yyzzz_yyyy[k] = -g_z_0_yzzz_yyyy[k] * cd_y[k] + g_z_0_yzzz_yyyyy[k];

                g_z_0_yyzzz_yyyz[k] = -g_z_0_yzzz_yyyz[k] * cd_y[k] + g_z_0_yzzz_yyyyz[k];

                g_z_0_yyzzz_yyzz[k] = -g_z_0_yzzz_yyzz[k] * cd_y[k] + g_z_0_yzzz_yyyzz[k];

                g_z_0_yyzzz_yzzz[k] = -g_z_0_yzzz_yzzz[k] * cd_y[k] + g_z_0_yzzz_yyzzz[k];

                g_z_0_yyzzz_zzzz[k] = -g_z_0_yzzz_zzzz[k] * cd_y[k] + g_z_0_yzzz_yzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 285);

            auto g_z_0_yzzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 286);

            auto g_z_0_yzzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 287);

            auto g_z_0_yzzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 288);

            auto g_z_0_yzzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 289);

            auto g_z_0_yzzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 290);

            auto g_z_0_yzzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 291);

            auto g_z_0_yzzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 292);

            auto g_z_0_yzzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 293);

            auto g_z_0_yzzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 294);

            auto g_z_0_yzzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 295);

            auto g_z_0_yzzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 296);

            auto g_z_0_yzzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 297);

            auto g_z_0_yzzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 298);

            auto g_z_0_yzzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 299);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzz_xxxx, g_z_0_yzzzz_xxxy, g_z_0_yzzzz_xxxz, g_z_0_yzzzz_xxyy, g_z_0_yzzzz_xxyz, g_z_0_yzzzz_xxzz, g_z_0_yzzzz_xyyy, g_z_0_yzzzz_xyyz, g_z_0_yzzzz_xyzz, g_z_0_yzzzz_xzzz, g_z_0_yzzzz_yyyy, g_z_0_yzzzz_yyyz, g_z_0_yzzzz_yyzz, g_z_0_yzzzz_yzzz, g_z_0_yzzzz_zzzz, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyyy, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_xxxx[k] = -g_z_0_zzzz_xxxx[k] * cd_y[k] + g_z_0_zzzz_xxxxy[k];

                g_z_0_yzzzz_xxxy[k] = -g_z_0_zzzz_xxxy[k] * cd_y[k] + g_z_0_zzzz_xxxyy[k];

                g_z_0_yzzzz_xxxz[k] = -g_z_0_zzzz_xxxz[k] * cd_y[k] + g_z_0_zzzz_xxxyz[k];

                g_z_0_yzzzz_xxyy[k] = -g_z_0_zzzz_xxyy[k] * cd_y[k] + g_z_0_zzzz_xxyyy[k];

                g_z_0_yzzzz_xxyz[k] = -g_z_0_zzzz_xxyz[k] * cd_y[k] + g_z_0_zzzz_xxyyz[k];

                g_z_0_yzzzz_xxzz[k] = -g_z_0_zzzz_xxzz[k] * cd_y[k] + g_z_0_zzzz_xxyzz[k];

                g_z_0_yzzzz_xyyy[k] = -g_z_0_zzzz_xyyy[k] * cd_y[k] + g_z_0_zzzz_xyyyy[k];

                g_z_0_yzzzz_xyyz[k] = -g_z_0_zzzz_xyyz[k] * cd_y[k] + g_z_0_zzzz_xyyyz[k];

                g_z_0_yzzzz_xyzz[k] = -g_z_0_zzzz_xyzz[k] * cd_y[k] + g_z_0_zzzz_xyyzz[k];

                g_z_0_yzzzz_xzzz[k] = -g_z_0_zzzz_xzzz[k] * cd_y[k] + g_z_0_zzzz_xyzzz[k];

                g_z_0_yzzzz_yyyy[k] = -g_z_0_zzzz_yyyy[k] * cd_y[k] + g_z_0_zzzz_yyyyy[k];

                g_z_0_yzzzz_yyyz[k] = -g_z_0_zzzz_yyyz[k] * cd_y[k] + g_z_0_zzzz_yyyyz[k];

                g_z_0_yzzzz_yyzz[k] = -g_z_0_zzzz_yyzz[k] * cd_y[k] + g_z_0_zzzz_yyyzz[k];

                g_z_0_yzzzz_yzzz[k] = -g_z_0_zzzz_yzzz[k] * cd_y[k] + g_z_0_zzzz_yyzzz[k];

                g_z_0_yzzzz_zzzz[k] = -g_z_0_zzzz_zzzz[k] * cd_y[k] + g_z_0_zzzz_yzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 300);

            auto g_z_0_zzzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 301);

            auto g_z_0_zzzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 302);

            auto g_z_0_zzzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 303);

            auto g_z_0_zzzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 304);

            auto g_z_0_zzzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 305);

            auto g_z_0_zzzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 306);

            auto g_z_0_zzzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 307);

            auto g_z_0_zzzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 308);

            auto g_z_0_zzzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 309);

            auto g_z_0_zzzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 310);

            auto g_z_0_zzzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 311);

            auto g_z_0_zzzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 312);

            auto g_z_0_zzzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 313);

            auto g_z_0_zzzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_z, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_zzzz, g_z_0_zzzz_zzzzz, g_z_0_zzzzz_xxxx, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_yyyy, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_zzzz, g_zzzz_xxxx, g_zzzz_xxxy, g_zzzz_xxxz, g_zzzz_xxyy, g_zzzz_xxyz, g_zzzz_xxzz, g_zzzz_xyyy, g_zzzz_xyyz, g_zzzz_xyzz, g_zzzz_xzzz, g_zzzz_yyyy, g_zzzz_yyyz, g_zzzz_yyzz, g_zzzz_yzzz, g_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_xxxx[k] = -g_zzzz_xxxx[k] - g_z_0_zzzz_xxxx[k] * cd_z[k] + g_z_0_zzzz_xxxxz[k];

                g_z_0_zzzzz_xxxy[k] = -g_zzzz_xxxy[k] - g_z_0_zzzz_xxxy[k] * cd_z[k] + g_z_0_zzzz_xxxyz[k];

                g_z_0_zzzzz_xxxz[k] = -g_zzzz_xxxz[k] - g_z_0_zzzz_xxxz[k] * cd_z[k] + g_z_0_zzzz_xxxzz[k];

                g_z_0_zzzzz_xxyy[k] = -g_zzzz_xxyy[k] - g_z_0_zzzz_xxyy[k] * cd_z[k] + g_z_0_zzzz_xxyyz[k];

                g_z_0_zzzzz_xxyz[k] = -g_zzzz_xxyz[k] - g_z_0_zzzz_xxyz[k] * cd_z[k] + g_z_0_zzzz_xxyzz[k];

                g_z_0_zzzzz_xxzz[k] = -g_zzzz_xxzz[k] - g_z_0_zzzz_xxzz[k] * cd_z[k] + g_z_0_zzzz_xxzzz[k];

                g_z_0_zzzzz_xyyy[k] = -g_zzzz_xyyy[k] - g_z_0_zzzz_xyyy[k] * cd_z[k] + g_z_0_zzzz_xyyyz[k];

                g_z_0_zzzzz_xyyz[k] = -g_zzzz_xyyz[k] - g_z_0_zzzz_xyyz[k] * cd_z[k] + g_z_0_zzzz_xyyzz[k];

                g_z_0_zzzzz_xyzz[k] = -g_zzzz_xyzz[k] - g_z_0_zzzz_xyzz[k] * cd_z[k] + g_z_0_zzzz_xyzzz[k];

                g_z_0_zzzzz_xzzz[k] = -g_zzzz_xzzz[k] - g_z_0_zzzz_xzzz[k] * cd_z[k] + g_z_0_zzzz_xzzzz[k];

                g_z_0_zzzzz_yyyy[k] = -g_zzzz_yyyy[k] - g_z_0_zzzz_yyyy[k] * cd_z[k] + g_z_0_zzzz_yyyyz[k];

                g_z_0_zzzzz_yyyz[k] = -g_zzzz_yyyz[k] - g_z_0_zzzz_yyyz[k] * cd_z[k] + g_z_0_zzzz_yyyzz[k];

                g_z_0_zzzzz_yyzz[k] = -g_zzzz_yyzz[k] - g_z_0_zzzz_yyzz[k] * cd_z[k] + g_z_0_zzzz_yyzzz[k];

                g_z_0_zzzzz_yzzz[k] = -g_zzzz_yzzz[k] - g_z_0_zzzz_yzzz[k] * cd_z[k] + g_z_0_zzzz_yzzzz[k];

                g_z_0_zzzzz_zzzz[k] = -g_zzzz_zzzz[k] - g_z_0_zzzz_zzzz[k] * cd_z[k] + g_z_0_zzzz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

