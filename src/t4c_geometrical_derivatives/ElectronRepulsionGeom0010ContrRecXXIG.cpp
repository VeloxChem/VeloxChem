#include "ElectronRepulsionGeom0010ContrRecXXIG.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxig(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxig,
                                            const size_t idx_xxhg,
                                            const size_t idx_geom_10_xxhg,
                                            const size_t idx_geom_10_xxhh,
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
            /// Set up components of auxilary buffer : SSHG

            const auto hg_off = idx_xxhg + (i * bcomps + j) * 315;

            auto g_xxxxx_xxxx = cbuffer.data(hg_off + 0);

            auto g_xxxxx_xxxy = cbuffer.data(hg_off + 1);

            auto g_xxxxx_xxxz = cbuffer.data(hg_off + 2);

            auto g_xxxxx_xxyy = cbuffer.data(hg_off + 3);

            auto g_xxxxx_xxyz = cbuffer.data(hg_off + 4);

            auto g_xxxxx_xxzz = cbuffer.data(hg_off + 5);

            auto g_xxxxx_xyyy = cbuffer.data(hg_off + 6);

            auto g_xxxxx_xyyz = cbuffer.data(hg_off + 7);

            auto g_xxxxx_xyzz = cbuffer.data(hg_off + 8);

            auto g_xxxxx_xzzz = cbuffer.data(hg_off + 9);

            auto g_xxxxx_yyyy = cbuffer.data(hg_off + 10);

            auto g_xxxxx_yyyz = cbuffer.data(hg_off + 11);

            auto g_xxxxx_yyzz = cbuffer.data(hg_off + 12);

            auto g_xxxxx_yzzz = cbuffer.data(hg_off + 13);

            auto g_xxxxx_zzzz = cbuffer.data(hg_off + 14);

            auto g_yyyyy_xxxx = cbuffer.data(hg_off + 225);

            auto g_yyyyy_xxxy = cbuffer.data(hg_off + 226);

            auto g_yyyyy_xxxz = cbuffer.data(hg_off + 227);

            auto g_yyyyy_xxyy = cbuffer.data(hg_off + 228);

            auto g_yyyyy_xxyz = cbuffer.data(hg_off + 229);

            auto g_yyyyy_xxzz = cbuffer.data(hg_off + 230);

            auto g_yyyyy_xyyy = cbuffer.data(hg_off + 231);

            auto g_yyyyy_xyyz = cbuffer.data(hg_off + 232);

            auto g_yyyyy_xyzz = cbuffer.data(hg_off + 233);

            auto g_yyyyy_xzzz = cbuffer.data(hg_off + 234);

            auto g_yyyyy_yyyy = cbuffer.data(hg_off + 235);

            auto g_yyyyy_yyyz = cbuffer.data(hg_off + 236);

            auto g_yyyyy_yyzz = cbuffer.data(hg_off + 237);

            auto g_yyyyy_yzzz = cbuffer.data(hg_off + 238);

            auto g_yyyyy_zzzz = cbuffer.data(hg_off + 239);

            auto g_zzzzz_xxxx = cbuffer.data(hg_off + 300);

            auto g_zzzzz_xxxy = cbuffer.data(hg_off + 301);

            auto g_zzzzz_xxxz = cbuffer.data(hg_off + 302);

            auto g_zzzzz_xxyy = cbuffer.data(hg_off + 303);

            auto g_zzzzz_xxyz = cbuffer.data(hg_off + 304);

            auto g_zzzzz_xxzz = cbuffer.data(hg_off + 305);

            auto g_zzzzz_xyyy = cbuffer.data(hg_off + 306);

            auto g_zzzzz_xyyz = cbuffer.data(hg_off + 307);

            auto g_zzzzz_xyzz = cbuffer.data(hg_off + 308);

            auto g_zzzzz_xzzz = cbuffer.data(hg_off + 309);

            auto g_zzzzz_yyyy = cbuffer.data(hg_off + 310);

            auto g_zzzzz_yyyz = cbuffer.data(hg_off + 311);

            auto g_zzzzz_yyzz = cbuffer.data(hg_off + 312);

            auto g_zzzzz_yzzz = cbuffer.data(hg_off + 313);

            auto g_zzzzz_zzzz = cbuffer.data(hg_off + 314);

            /// Set up components of auxilary buffer : SSHG

            const auto hg_geom_10_off = idx_geom_10_xxhg + (i * bcomps + j) * 315;

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

            /// Set up components of auxilary buffer : SSHH

            const auto hh_geom_10_off = idx_geom_10_xxhh + (i * bcomps + j) * 441;

            auto g_x_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_x_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 329);

            auto g_x_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 335);

            auto g_x_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 344);

            auto g_x_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 356);

            auto g_x_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 359);

            auto g_x_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 360);

            auto g_x_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 361);

            auto g_x_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 362);

            auto g_x_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 363);

            auto g_x_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 364);

            auto g_x_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 365);

            auto g_x_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 366);

            auto g_x_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 367);

            auto g_x_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 368);

            auto g_x_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 369);

            auto g_x_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 370);

            auto g_x_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 371);

            auto g_x_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 372);

            auto g_x_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 373);

            auto g_x_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 374);

            auto g_x_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 375);

            auto g_x_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 376);

            auto g_x_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 377);

            auto g_x_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 378);

            auto g_x_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 379);

            auto g_x_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 380);

            auto g_x_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 381);

            auto g_x_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 382);

            auto g_x_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 383);

            auto g_x_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 384);

            auto g_x_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 385);

            auto g_x_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 386);

            auto g_x_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 387);

            auto g_x_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 388);

            auto g_x_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 389);

            auto g_x_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 390);

            auto g_x_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 391);

            auto g_x_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 392);

            auto g_x_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 393);

            auto g_x_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 394);

            auto g_x_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 395);

            auto g_x_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 396);

            auto g_x_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 397);

            auto g_x_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 398);

            auto g_x_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 399);

            auto g_x_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 400);

            auto g_x_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 401);

            auto g_x_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 402);

            auto g_x_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 403);

            auto g_x_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 404);

            auto g_x_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 405);

            auto g_x_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 406);

            auto g_x_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 407);

            auto g_x_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 408);

            auto g_x_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 409);

            auto g_x_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 410);

            auto g_x_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 411);

            auto g_x_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 412);

            auto g_x_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 413);

            auto g_x_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 414);

            auto g_x_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 415);

            auto g_x_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 416);

            auto g_x_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 417);

            auto g_x_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 418);

            auto g_x_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 419);

            auto g_x_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 420);

            auto g_x_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 421);

            auto g_x_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 422);

            auto g_x_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 423);

            auto g_x_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 424);

            auto g_x_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 425);

            auto g_x_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 426);

            auto g_x_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 427);

            auto g_x_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 428);

            auto g_x_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 429);

            auto g_x_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 430);

            auto g_x_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 431);

            auto g_x_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 432);

            auto g_x_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 433);

            auto g_x_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 434);

            auto g_x_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 435);

            auto g_x_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 436);

            auto g_x_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 437);

            auto g_x_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 438);

            auto g_x_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 439);

            auto g_x_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 440);

            auto g_y_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 5);

            auto g_y_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 6);

            auto g_y_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 7);

            auto g_y_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 8);

            auto g_y_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 9);

            auto g_y_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 10);

            auto g_y_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 11);

            auto g_y_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 12);

            auto g_y_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 13);

            auto g_y_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 14);

            auto g_y_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 15);

            auto g_y_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 16);

            auto g_y_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 17);

            auto g_y_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 18);

            auto g_y_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 19);

            auto g_y_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 20);

            auto g_y_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 21);

            auto g_y_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 22);

            auto g_y_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 23);

            auto g_y_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 24);

            auto g_y_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 25);

            auto g_y_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 26);

            auto g_y_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 27);

            auto g_y_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 28);

            auto g_y_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 29);

            auto g_y_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 30);

            auto g_y_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 31);

            auto g_y_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 32);

            auto g_y_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 33);

            auto g_y_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 34);

            auto g_y_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 35);

            auto g_y_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 36);

            auto g_y_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 37);

            auto g_y_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 38);

            auto g_y_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 39);

            auto g_y_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 40);

            auto g_y_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 41);

            auto g_y_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 42);

            auto g_y_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 43);

            auto g_y_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 44);

            auto g_y_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 45);

            auto g_y_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 46);

            auto g_y_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 47);

            auto g_y_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 48);

            auto g_y_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 49);

            auto g_y_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 50);

            auto g_y_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 51);

            auto g_y_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 52);

            auto g_y_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 53);

            auto g_y_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 54);

            auto g_y_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 55);

            auto g_y_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 56);

            auto g_y_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 57);

            auto g_y_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 58);

            auto g_y_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 59);

            auto g_y_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 60);

            auto g_y_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 61);

            auto g_y_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 62);

            auto g_y_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 63);

            auto g_y_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 64);

            auto g_y_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 65);

            auto g_y_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 66);

            auto g_y_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 67);

            auto g_y_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 68);

            auto g_y_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 69);

            auto g_y_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 70);

            auto g_y_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 71);

            auto g_y_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 72);

            auto g_y_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 73);

            auto g_y_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 74);

            auto g_y_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 75);

            auto g_y_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 76);

            auto g_y_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 77);

            auto g_y_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 78);

            auto g_y_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 79);

            auto g_y_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 80);

            auto g_y_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 81);

            auto g_y_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 82);

            auto g_y_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 83);

            auto g_y_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 84);

            auto g_y_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 85);

            auto g_y_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 86);

            auto g_y_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 87);

            auto g_y_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 88);

            auto g_y_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 89);

            auto g_y_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 90);

            auto g_y_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 91);

            auto g_y_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 92);

            auto g_y_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 93);

            auto g_y_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 94);

            auto g_y_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 95);

            auto g_y_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 96);

            auto g_y_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 97);

            auto g_y_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 98);

            auto g_y_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 99);

            auto g_y_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 100);

            auto g_y_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 101);

            auto g_y_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 102);

            auto g_y_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 103);

            auto g_y_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 104);

            auto g_y_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 105);

            auto g_y_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 106);

            auto g_y_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 107);

            auto g_y_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 108);

            auto g_y_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 109);

            auto g_y_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 110);

            auto g_y_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 111);

            auto g_y_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 112);

            auto g_y_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 113);

            auto g_y_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 114);

            auto g_y_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 115);

            auto g_y_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 116);

            auto g_y_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 117);

            auto g_y_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 118);

            auto g_y_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 119);

            auto g_y_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 120);

            auto g_y_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 121);

            auto g_y_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 122);

            auto g_y_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 123);

            auto g_y_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 124);

            auto g_y_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 125);

            auto g_y_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 126);

            auto g_y_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 127);

            auto g_y_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 128);

            auto g_y_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 129);

            auto g_y_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 130);

            auto g_y_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 131);

            auto g_y_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 132);

            auto g_y_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 133);

            auto g_y_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 134);

            auto g_y_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 135);

            auto g_y_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 136);

            auto g_y_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 137);

            auto g_y_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 138);

            auto g_y_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 139);

            auto g_y_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 140);

            auto g_y_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 141);

            auto g_y_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 142);

            auto g_y_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 143);

            auto g_y_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 144);

            auto g_y_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 145);

            auto g_y_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 146);

            auto g_y_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 147);

            auto g_y_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 148);

            auto g_y_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 149);

            auto g_y_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 150);

            auto g_y_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 151);

            auto g_y_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 152);

            auto g_y_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 153);

            auto g_y_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 154);

            auto g_y_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 155);

            auto g_y_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 156);

            auto g_y_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 157);

            auto g_y_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 158);

            auto g_y_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 159);

            auto g_y_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 160);

            auto g_y_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 161);

            auto g_y_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 162);

            auto g_y_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 163);

            auto g_y_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 164);

            auto g_y_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 165);

            auto g_y_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 166);

            auto g_y_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 167);

            auto g_y_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 168);

            auto g_y_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 169);

            auto g_y_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 170);

            auto g_y_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 171);

            auto g_y_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 172);

            auto g_y_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 173);

            auto g_y_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 174);

            auto g_y_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 175);

            auto g_y_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 176);

            auto g_y_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 177);

            auto g_y_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 178);

            auto g_y_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 179);

            auto g_y_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 180);

            auto g_y_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 181);

            auto g_y_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 182);

            auto g_y_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 183);

            auto g_y_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 184);

            auto g_y_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 185);

            auto g_y_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 186);

            auto g_y_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 187);

            auto g_y_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 188);

            auto g_y_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 189);

            auto g_y_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 190);

            auto g_y_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 191);

            auto g_y_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 192);

            auto g_y_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 193);

            auto g_y_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 194);

            auto g_y_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 195);

            auto g_y_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 196);

            auto g_y_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 197);

            auto g_y_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 198);

            auto g_y_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 199);

            auto g_y_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 200);

            auto g_y_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 201);

            auto g_y_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 202);

            auto g_y_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 203);

            auto g_y_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 204);

            auto g_y_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 205);

            auto g_y_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 206);

            auto g_y_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 207);

            auto g_y_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 208);

            auto g_y_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 209);

            auto g_y_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 210);

            auto g_y_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 211);

            auto g_y_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 212);

            auto g_y_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 213);

            auto g_y_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 214);

            auto g_y_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 215);

            auto g_y_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 216);

            auto g_y_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 217);

            auto g_y_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 218);

            auto g_y_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 219);

            auto g_y_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 220);

            auto g_y_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 221);

            auto g_y_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 222);

            auto g_y_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 223);

            auto g_y_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 224);

            auto g_y_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 225);

            auto g_y_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 226);

            auto g_y_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 227);

            auto g_y_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 228);

            auto g_y_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 229);

            auto g_y_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 230);

            auto g_y_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 231);

            auto g_y_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 232);

            auto g_y_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 233);

            auto g_y_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 234);

            auto g_y_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 235);

            auto g_y_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 236);

            auto g_y_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 237);

            auto g_y_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 238);

            auto g_y_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 239);

            auto g_y_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 240);

            auto g_y_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 241);

            auto g_y_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 242);

            auto g_y_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 243);

            auto g_y_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 244);

            auto g_y_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 245);

            auto g_y_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 246);

            auto g_y_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 247);

            auto g_y_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 248);

            auto g_y_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 249);

            auto g_y_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 250);

            auto g_y_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 251);

            auto g_y_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 252);

            auto g_y_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 253);

            auto g_y_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 254);

            auto g_y_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 255);

            auto g_y_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 256);

            auto g_y_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 257);

            auto g_y_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 258);

            auto g_y_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 259);

            auto g_y_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 260);

            auto g_y_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 261);

            auto g_y_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 262);

            auto g_y_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 263);

            auto g_y_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 264);

            auto g_y_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 265);

            auto g_y_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 266);

            auto g_y_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 267);

            auto g_y_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 268);

            auto g_y_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 269);

            auto g_y_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 270);

            auto g_y_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 271);

            auto g_y_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 272);

            auto g_y_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 273);

            auto g_y_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 274);

            auto g_y_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 275);

            auto g_y_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 276);

            auto g_y_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 277);

            auto g_y_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 278);

            auto g_y_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 279);

            auto g_y_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 280);

            auto g_y_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 281);

            auto g_y_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 282);

            auto g_y_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 283);

            auto g_y_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 284);

            auto g_y_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 285);

            auto g_y_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 286);

            auto g_y_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 287);

            auto g_y_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 288);

            auto g_y_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 289);

            auto g_y_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 290);

            auto g_y_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 291);

            auto g_y_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 292);

            auto g_y_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 293);

            auto g_y_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 294);

            auto g_y_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 295);

            auto g_y_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 296);

            auto g_y_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 297);

            auto g_y_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 298);

            auto g_y_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 299);

            auto g_y_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 300);

            auto g_y_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 301);

            auto g_y_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 302);

            auto g_y_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 303);

            auto g_y_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 304);

            auto g_y_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 305);

            auto g_y_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 306);

            auto g_y_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 307);

            auto g_y_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 308);

            auto g_y_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 309);

            auto g_y_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 310);

            auto g_y_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 311);

            auto g_y_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 312);

            auto g_y_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 313);

            auto g_y_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 314);

            auto g_y_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 315);

            auto g_y_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 316);

            auto g_y_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 317);

            auto g_y_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 318);

            auto g_y_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 319);

            auto g_y_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 320);

            auto g_y_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 321);

            auto g_y_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 322);

            auto g_y_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 323);

            auto g_y_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 324);

            auto g_y_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 325);

            auto g_y_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 326);

            auto g_y_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 327);

            auto g_y_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 328);

            auto g_y_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 329);

            auto g_y_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 330);

            auto g_y_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 331);

            auto g_y_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 332);

            auto g_y_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 333);

            auto g_y_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 334);

            auto g_y_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 335);

            auto g_y_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 336);

            auto g_y_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 337);

            auto g_y_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 338);

            auto g_y_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 339);

            auto g_y_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 340);

            auto g_y_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 341);

            auto g_y_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 342);

            auto g_y_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 343);

            auto g_y_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 344);

            auto g_y_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 345);

            auto g_y_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 346);

            auto g_y_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 347);

            auto g_y_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 348);

            auto g_y_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 349);

            auto g_y_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 350);

            auto g_y_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 351);

            auto g_y_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 352);

            auto g_y_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 353);

            auto g_y_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 354);

            auto g_y_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 355);

            auto g_y_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 356);

            auto g_y_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 357);

            auto g_y_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 358);

            auto g_y_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 359);

            auto g_y_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 360);

            auto g_y_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 361);

            auto g_y_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 362);

            auto g_y_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 363);

            auto g_y_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 364);

            auto g_y_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 365);

            auto g_y_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 366);

            auto g_y_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 367);

            auto g_y_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 368);

            auto g_y_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 369);

            auto g_y_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 370);

            auto g_y_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 371);

            auto g_y_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 372);

            auto g_y_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 373);

            auto g_y_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 374);

            auto g_y_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 375);

            auto g_y_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 376);

            auto g_y_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 377);

            auto g_y_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 378);

            auto g_y_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 379);

            auto g_y_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 380);

            auto g_y_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 381);

            auto g_y_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 382);

            auto g_y_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 383);

            auto g_y_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 384);

            auto g_y_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 385);

            auto g_y_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 386);

            auto g_y_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 387);

            auto g_y_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 388);

            auto g_y_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 389);

            auto g_y_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 390);

            auto g_y_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 391);

            auto g_y_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 392);

            auto g_y_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 393);

            auto g_y_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 394);

            auto g_y_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 395);

            auto g_y_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 396);

            auto g_y_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 397);

            auto g_y_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 398);

            auto g_y_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 399);

            auto g_y_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 400);

            auto g_y_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 401);

            auto g_y_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 402);

            auto g_y_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 403);

            auto g_y_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 404);

            auto g_y_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 405);

            auto g_y_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 406);

            auto g_y_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 407);

            auto g_y_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 408);

            auto g_y_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 409);

            auto g_y_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 410);

            auto g_y_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 411);

            auto g_y_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 412);

            auto g_y_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 413);

            auto g_y_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 414);

            auto g_y_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 415);

            auto g_y_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 416);

            auto g_y_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 417);

            auto g_y_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 418);

            auto g_y_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 419);

            auto g_y_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 420);

            auto g_y_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 421);

            auto g_y_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 422);

            auto g_y_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 423);

            auto g_y_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 424);

            auto g_y_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 425);

            auto g_y_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 426);

            auto g_y_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 427);

            auto g_y_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 428);

            auto g_y_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 429);

            auto g_y_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 430);

            auto g_y_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 431);

            auto g_y_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 432);

            auto g_y_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 433);

            auto g_y_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 434);

            auto g_y_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 435);

            auto g_y_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 436);

            auto g_y_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 437);

            auto g_y_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 438);

            auto g_y_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 439);

            auto g_y_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 440);

            auto g_z_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 5);

            auto g_z_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 6);

            auto g_z_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 7);

            auto g_z_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 8);

            auto g_z_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 9);

            auto g_z_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 10);

            auto g_z_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 11);

            auto g_z_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 12);

            auto g_z_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 13);

            auto g_z_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 14);

            auto g_z_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 15);

            auto g_z_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 16);

            auto g_z_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 17);

            auto g_z_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 18);

            auto g_z_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 19);

            auto g_z_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 20);

            auto g_z_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 21);

            auto g_z_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 22);

            auto g_z_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 23);

            auto g_z_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 24);

            auto g_z_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 25);

            auto g_z_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 26);

            auto g_z_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 27);

            auto g_z_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 28);

            auto g_z_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 29);

            auto g_z_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 30);

            auto g_z_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 31);

            auto g_z_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 32);

            auto g_z_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 33);

            auto g_z_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 34);

            auto g_z_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 35);

            auto g_z_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 36);

            auto g_z_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 37);

            auto g_z_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 38);

            auto g_z_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 39);

            auto g_z_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 40);

            auto g_z_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 41);

            auto g_z_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 42);

            auto g_z_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 43);

            auto g_z_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 44);

            auto g_z_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 45);

            auto g_z_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 46);

            auto g_z_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 47);

            auto g_z_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 48);

            auto g_z_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 49);

            auto g_z_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 50);

            auto g_z_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 51);

            auto g_z_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 52);

            auto g_z_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 53);

            auto g_z_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 54);

            auto g_z_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 55);

            auto g_z_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 56);

            auto g_z_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 57);

            auto g_z_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 58);

            auto g_z_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 59);

            auto g_z_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 60);

            auto g_z_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 61);

            auto g_z_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 62);

            auto g_z_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 63);

            auto g_z_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 64);

            auto g_z_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 65);

            auto g_z_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 66);

            auto g_z_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 67);

            auto g_z_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 68);

            auto g_z_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 69);

            auto g_z_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 70);

            auto g_z_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 71);

            auto g_z_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 72);

            auto g_z_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 73);

            auto g_z_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 74);

            auto g_z_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 75);

            auto g_z_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 76);

            auto g_z_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 77);

            auto g_z_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 78);

            auto g_z_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 79);

            auto g_z_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 80);

            auto g_z_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 81);

            auto g_z_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 82);

            auto g_z_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 83);

            auto g_z_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 84);

            auto g_z_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 85);

            auto g_z_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 86);

            auto g_z_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 87);

            auto g_z_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 88);

            auto g_z_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 89);

            auto g_z_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 90);

            auto g_z_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 91);

            auto g_z_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 92);

            auto g_z_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 93);

            auto g_z_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 94);

            auto g_z_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 95);

            auto g_z_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 96);

            auto g_z_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 97);

            auto g_z_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 98);

            auto g_z_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 99);

            auto g_z_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 100);

            auto g_z_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 101);

            auto g_z_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 102);

            auto g_z_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 103);

            auto g_z_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 104);

            auto g_z_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 105);

            auto g_z_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 106);

            auto g_z_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 107);

            auto g_z_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 108);

            auto g_z_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 109);

            auto g_z_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 110);

            auto g_z_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 111);

            auto g_z_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 112);

            auto g_z_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 113);

            auto g_z_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 114);

            auto g_z_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 115);

            auto g_z_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 116);

            auto g_z_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 117);

            auto g_z_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 118);

            auto g_z_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 119);

            auto g_z_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 120);

            auto g_z_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 121);

            auto g_z_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 122);

            auto g_z_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 123);

            auto g_z_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 124);

            auto g_z_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 125);

            auto g_z_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 126);

            auto g_z_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 127);

            auto g_z_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 128);

            auto g_z_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 129);

            auto g_z_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 130);

            auto g_z_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 131);

            auto g_z_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 132);

            auto g_z_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 133);

            auto g_z_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 134);

            auto g_z_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 135);

            auto g_z_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 136);

            auto g_z_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 137);

            auto g_z_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 138);

            auto g_z_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 139);

            auto g_z_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 140);

            auto g_z_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 141);

            auto g_z_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 142);

            auto g_z_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 143);

            auto g_z_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 144);

            auto g_z_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 145);

            auto g_z_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 146);

            auto g_z_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 147);

            auto g_z_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 148);

            auto g_z_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 149);

            auto g_z_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 150);

            auto g_z_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 151);

            auto g_z_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 152);

            auto g_z_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 153);

            auto g_z_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 154);

            auto g_z_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 155);

            auto g_z_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 156);

            auto g_z_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 157);

            auto g_z_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 158);

            auto g_z_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 159);

            auto g_z_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 160);

            auto g_z_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 161);

            auto g_z_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 162);

            auto g_z_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 163);

            auto g_z_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 164);

            auto g_z_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 165);

            auto g_z_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 166);

            auto g_z_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 167);

            auto g_z_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 168);

            auto g_z_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 169);

            auto g_z_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 170);

            auto g_z_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 171);

            auto g_z_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 172);

            auto g_z_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 173);

            auto g_z_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 174);

            auto g_z_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 175);

            auto g_z_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 176);

            auto g_z_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 177);

            auto g_z_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 178);

            auto g_z_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 179);

            auto g_z_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 180);

            auto g_z_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 181);

            auto g_z_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 182);

            auto g_z_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 183);

            auto g_z_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 184);

            auto g_z_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 185);

            auto g_z_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 186);

            auto g_z_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 187);

            auto g_z_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 188);

            auto g_z_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 189);

            auto g_z_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 190);

            auto g_z_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 191);

            auto g_z_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 192);

            auto g_z_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 193);

            auto g_z_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 194);

            auto g_z_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 195);

            auto g_z_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 196);

            auto g_z_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 197);

            auto g_z_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 198);

            auto g_z_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 199);

            auto g_z_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 200);

            auto g_z_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 201);

            auto g_z_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 202);

            auto g_z_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 203);

            auto g_z_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 204);

            auto g_z_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 205);

            auto g_z_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 206);

            auto g_z_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 207);

            auto g_z_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 208);

            auto g_z_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 209);

            auto g_z_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 210);

            auto g_z_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 211);

            auto g_z_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 212);

            auto g_z_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 213);

            auto g_z_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 214);

            auto g_z_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 215);

            auto g_z_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 216);

            auto g_z_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 217);

            auto g_z_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 218);

            auto g_z_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 219);

            auto g_z_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 220);

            auto g_z_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 221);

            auto g_z_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 222);

            auto g_z_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 223);

            auto g_z_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 224);

            auto g_z_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 225);

            auto g_z_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 226);

            auto g_z_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 227);

            auto g_z_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 228);

            auto g_z_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 229);

            auto g_z_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 230);

            auto g_z_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 231);

            auto g_z_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 232);

            auto g_z_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 233);

            auto g_z_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 234);

            auto g_z_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 235);

            auto g_z_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 236);

            auto g_z_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 237);

            auto g_z_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 238);

            auto g_z_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 239);

            auto g_z_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 240);

            auto g_z_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 241);

            auto g_z_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 242);

            auto g_z_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 243);

            auto g_z_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 244);

            auto g_z_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 245);

            auto g_z_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 246);

            auto g_z_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 247);

            auto g_z_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 248);

            auto g_z_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 249);

            auto g_z_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 250);

            auto g_z_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 251);

            auto g_z_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 252);

            auto g_z_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 253);

            auto g_z_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 254);

            auto g_z_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 255);

            auto g_z_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 256);

            auto g_z_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 257);

            auto g_z_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 258);

            auto g_z_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 259);

            auto g_z_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 260);

            auto g_z_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 261);

            auto g_z_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 262);

            auto g_z_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 263);

            auto g_z_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 264);

            auto g_z_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 265);

            auto g_z_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 266);

            auto g_z_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 267);

            auto g_z_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 268);

            auto g_z_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 269);

            auto g_z_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 270);

            auto g_z_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 271);

            auto g_z_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 272);

            auto g_z_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 273);

            auto g_z_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 274);

            auto g_z_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 275);

            auto g_z_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 276);

            auto g_z_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 277);

            auto g_z_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 278);

            auto g_z_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 279);

            auto g_z_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 280);

            auto g_z_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 281);

            auto g_z_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 282);

            auto g_z_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 283);

            auto g_z_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 284);

            auto g_z_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 285);

            auto g_z_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 286);

            auto g_z_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 287);

            auto g_z_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 288);

            auto g_z_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 289);

            auto g_z_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 290);

            auto g_z_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 291);

            auto g_z_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 292);

            auto g_z_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 293);

            auto g_z_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 294);

            auto g_z_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 295);

            auto g_z_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 296);

            auto g_z_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 297);

            auto g_z_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 298);

            auto g_z_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 299);

            auto g_z_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 300);

            auto g_z_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 301);

            auto g_z_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 302);

            auto g_z_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 303);

            auto g_z_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 304);

            auto g_z_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 305);

            auto g_z_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 306);

            auto g_z_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 307);

            auto g_z_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 308);

            auto g_z_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 309);

            auto g_z_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 310);

            auto g_z_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 311);

            auto g_z_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 312);

            auto g_z_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 313);

            auto g_z_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 314);

            auto g_z_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 315);

            auto g_z_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 316);

            auto g_z_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 317);

            auto g_z_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 318);

            auto g_z_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 319);

            auto g_z_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 320);

            auto g_z_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 321);

            auto g_z_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 322);

            auto g_z_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 323);

            auto g_z_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 324);

            auto g_z_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 325);

            auto g_z_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 326);

            auto g_z_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 327);

            auto g_z_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 328);

            auto g_z_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 329);

            auto g_z_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 330);

            auto g_z_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 331);

            auto g_z_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 332);

            auto g_z_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 333);

            auto g_z_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 334);

            auto g_z_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 335);

            auto g_z_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 336);

            auto g_z_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 337);

            auto g_z_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 338);

            auto g_z_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 339);

            auto g_z_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 340);

            auto g_z_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 341);

            auto g_z_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 342);

            auto g_z_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 343);

            auto g_z_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 344);

            auto g_z_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 345);

            auto g_z_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 346);

            auto g_z_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 347);

            auto g_z_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 348);

            auto g_z_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 349);

            auto g_z_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 350);

            auto g_z_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 351);

            auto g_z_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 352);

            auto g_z_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 353);

            auto g_z_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 354);

            auto g_z_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 355);

            auto g_z_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 356);

            auto g_z_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 357);

            auto g_z_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 358);

            auto g_z_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 359);

            auto g_z_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 360);

            auto g_z_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 361);

            auto g_z_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 362);

            auto g_z_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 363);

            auto g_z_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 364);

            auto g_z_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 365);

            auto g_z_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 366);

            auto g_z_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 367);

            auto g_z_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 368);

            auto g_z_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 369);

            auto g_z_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 370);

            auto g_z_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 371);

            auto g_z_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 372);

            auto g_z_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 373);

            auto g_z_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 374);

            auto g_z_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 375);

            auto g_z_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 376);

            auto g_z_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 377);

            auto g_z_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 378);

            auto g_z_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 379);

            auto g_z_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 380);

            auto g_z_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 381);

            auto g_z_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 382);

            auto g_z_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 383);

            auto g_z_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 384);

            auto g_z_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 385);

            auto g_z_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 386);

            auto g_z_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 387);

            auto g_z_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 388);

            auto g_z_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 389);

            auto g_z_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 390);

            auto g_z_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 391);

            auto g_z_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 392);

            auto g_z_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 393);

            auto g_z_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 394);

            auto g_z_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 395);

            auto g_z_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 396);

            auto g_z_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 397);

            auto g_z_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 398);

            auto g_z_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 399);

            auto g_z_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 400);

            auto g_z_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 401);

            auto g_z_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 402);

            auto g_z_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 403);

            auto g_z_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 404);

            auto g_z_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 405);

            auto g_z_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 406);

            auto g_z_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 407);

            auto g_z_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 408);

            auto g_z_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 409);

            auto g_z_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 410);

            auto g_z_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 411);

            auto g_z_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 412);

            auto g_z_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 413);

            auto g_z_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 414);

            auto g_z_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 415);

            auto g_z_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 416);

            auto g_z_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 417);

            auto g_z_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 418);

            auto g_z_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 419);

            auto g_z_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 420);

            auto g_z_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 421);

            auto g_z_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 422);

            auto g_z_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 423);

            auto g_z_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 424);

            auto g_z_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 425);

            auto g_z_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 426);

            auto g_z_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 427);

            auto g_z_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 428);

            auto g_z_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 429);

            auto g_z_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 430);

            auto g_z_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 431);

            auto g_z_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 432);

            auto g_z_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 433);

            auto g_z_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 434);

            auto g_z_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 435);

            auto g_z_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 436);

            auto g_z_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 437);

            auto g_z_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 438);

            auto g_z_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 439);

            auto g_z_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 440);

            /// set up bra offset for contr_buffer_xxig

            const auto ig_geom_10_off = idx_geom_10_xxig + (i * bcomps + j) * 420;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxxx_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxxx_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxxx_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxxx_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxxx_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxxx_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxxx_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxxx_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxxx_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxxx_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxxx_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxxx_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxxx_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxxx_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_x_0_xxxxx_xxxx, g_x_0_xxxxx_xxxxx, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_yyyy, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_zzzz, g_x_0_xxxxxx_xxxx, g_x_0_xxxxxx_xxxy, g_x_0_xxxxxx_xxxz, g_x_0_xxxxxx_xxyy, g_x_0_xxxxxx_xxyz, g_x_0_xxxxxx_xxzz, g_x_0_xxxxxx_xyyy, g_x_0_xxxxxx_xyyz, g_x_0_xxxxxx_xyzz, g_x_0_xxxxxx_xzzz, g_x_0_xxxxxx_yyyy, g_x_0_xxxxxx_yyyz, g_x_0_xxxxxx_yyzz, g_x_0_xxxxxx_yzzz, g_x_0_xxxxxx_zzzz, g_xxxxx_xxxx, g_xxxxx_xxxy, g_xxxxx_xxxz, g_xxxxx_xxyy, g_xxxxx_xxyz, g_xxxxx_xxzz, g_xxxxx_xyyy, g_xxxxx_xyyz, g_xxxxx_xyzz, g_xxxxx_xzzz, g_xxxxx_yyyy, g_xxxxx_yyyz, g_xxxxx_yyzz, g_xxxxx_yzzz, g_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_xxxx[k] = -g_xxxxx_xxxx[k] - g_x_0_xxxxx_xxxx[k] * cd_x[k] + g_x_0_xxxxx_xxxxx[k];

                g_x_0_xxxxxx_xxxy[k] = -g_xxxxx_xxxy[k] - g_x_0_xxxxx_xxxy[k] * cd_x[k] + g_x_0_xxxxx_xxxxy[k];

                g_x_0_xxxxxx_xxxz[k] = -g_xxxxx_xxxz[k] - g_x_0_xxxxx_xxxz[k] * cd_x[k] + g_x_0_xxxxx_xxxxz[k];

                g_x_0_xxxxxx_xxyy[k] = -g_xxxxx_xxyy[k] - g_x_0_xxxxx_xxyy[k] * cd_x[k] + g_x_0_xxxxx_xxxyy[k];

                g_x_0_xxxxxx_xxyz[k] = -g_xxxxx_xxyz[k] - g_x_0_xxxxx_xxyz[k] * cd_x[k] + g_x_0_xxxxx_xxxyz[k];

                g_x_0_xxxxxx_xxzz[k] = -g_xxxxx_xxzz[k] - g_x_0_xxxxx_xxzz[k] * cd_x[k] + g_x_0_xxxxx_xxxzz[k];

                g_x_0_xxxxxx_xyyy[k] = -g_xxxxx_xyyy[k] - g_x_0_xxxxx_xyyy[k] * cd_x[k] + g_x_0_xxxxx_xxyyy[k];

                g_x_0_xxxxxx_xyyz[k] = -g_xxxxx_xyyz[k] - g_x_0_xxxxx_xyyz[k] * cd_x[k] + g_x_0_xxxxx_xxyyz[k];

                g_x_0_xxxxxx_xyzz[k] = -g_xxxxx_xyzz[k] - g_x_0_xxxxx_xyzz[k] * cd_x[k] + g_x_0_xxxxx_xxyzz[k];

                g_x_0_xxxxxx_xzzz[k] = -g_xxxxx_xzzz[k] - g_x_0_xxxxx_xzzz[k] * cd_x[k] + g_x_0_xxxxx_xxzzz[k];

                g_x_0_xxxxxx_yyyy[k] = -g_xxxxx_yyyy[k] - g_x_0_xxxxx_yyyy[k] * cd_x[k] + g_x_0_xxxxx_xyyyy[k];

                g_x_0_xxxxxx_yyyz[k] = -g_xxxxx_yyyz[k] - g_x_0_xxxxx_yyyz[k] * cd_x[k] + g_x_0_xxxxx_xyyyz[k];

                g_x_0_xxxxxx_yyzz[k] = -g_xxxxx_yyzz[k] - g_x_0_xxxxx_yyzz[k] * cd_x[k] + g_x_0_xxxxx_xyyzz[k];

                g_x_0_xxxxxx_yzzz[k] = -g_xxxxx_yzzz[k] - g_x_0_xxxxx_yzzz[k] * cd_x[k] + g_x_0_xxxxx_xyzzz[k];

                g_x_0_xxxxxx_zzzz[k] = -g_xxxxx_zzzz[k] - g_x_0_xxxxx_zzzz[k] * cd_x[k] + g_x_0_xxxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxxy_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxxy_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxxxy_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxxy_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxxxy_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxxxy_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxxy_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxxy_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxxy_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxxy_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxxy_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxxy_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxxy_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxxy_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxx_xxxx, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_yyyy, g_x_0_xxxxx_yyyyy, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_zzzz, g_x_0_xxxxxy_xxxx, g_x_0_xxxxxy_xxxy, g_x_0_xxxxxy_xxxz, g_x_0_xxxxxy_xxyy, g_x_0_xxxxxy_xxyz, g_x_0_xxxxxy_xxzz, g_x_0_xxxxxy_xyyy, g_x_0_xxxxxy_xyyz, g_x_0_xxxxxy_xyzz, g_x_0_xxxxxy_xzzz, g_x_0_xxxxxy_yyyy, g_x_0_xxxxxy_yyyz, g_x_0_xxxxxy_yyzz, g_x_0_xxxxxy_yzzz, g_x_0_xxxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_xxxx[k] = -g_x_0_xxxxx_xxxx[k] * cd_y[k] + g_x_0_xxxxx_xxxxy[k];

                g_x_0_xxxxxy_xxxy[k] = -g_x_0_xxxxx_xxxy[k] * cd_y[k] + g_x_0_xxxxx_xxxyy[k];

                g_x_0_xxxxxy_xxxz[k] = -g_x_0_xxxxx_xxxz[k] * cd_y[k] + g_x_0_xxxxx_xxxyz[k];

                g_x_0_xxxxxy_xxyy[k] = -g_x_0_xxxxx_xxyy[k] * cd_y[k] + g_x_0_xxxxx_xxyyy[k];

                g_x_0_xxxxxy_xxyz[k] = -g_x_0_xxxxx_xxyz[k] * cd_y[k] + g_x_0_xxxxx_xxyyz[k];

                g_x_0_xxxxxy_xxzz[k] = -g_x_0_xxxxx_xxzz[k] * cd_y[k] + g_x_0_xxxxx_xxyzz[k];

                g_x_0_xxxxxy_xyyy[k] = -g_x_0_xxxxx_xyyy[k] * cd_y[k] + g_x_0_xxxxx_xyyyy[k];

                g_x_0_xxxxxy_xyyz[k] = -g_x_0_xxxxx_xyyz[k] * cd_y[k] + g_x_0_xxxxx_xyyyz[k];

                g_x_0_xxxxxy_xyzz[k] = -g_x_0_xxxxx_xyzz[k] * cd_y[k] + g_x_0_xxxxx_xyyzz[k];

                g_x_0_xxxxxy_xzzz[k] = -g_x_0_xxxxx_xzzz[k] * cd_y[k] + g_x_0_xxxxx_xyzzz[k];

                g_x_0_xxxxxy_yyyy[k] = -g_x_0_xxxxx_yyyy[k] * cd_y[k] + g_x_0_xxxxx_yyyyy[k];

                g_x_0_xxxxxy_yyyz[k] = -g_x_0_xxxxx_yyyz[k] * cd_y[k] + g_x_0_xxxxx_yyyyz[k];

                g_x_0_xxxxxy_yyzz[k] = -g_x_0_xxxxx_yyzz[k] * cd_y[k] + g_x_0_xxxxx_yyyzz[k];

                g_x_0_xxxxxy_yzzz[k] = -g_x_0_xxxxx_yzzz[k] * cd_y[k] + g_x_0_xxxxx_yyzzz[k];

                g_x_0_xxxxxy_zzzz[k] = -g_x_0_xxxxx_zzzz[k] * cd_y[k] + g_x_0_xxxxx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxxxz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxxxz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxxxz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxxxz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxxxz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxxxz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxxxz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxxxz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxxxz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxxxz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxxxz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxxxz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxxxz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxxxz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxx_xxxx, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_yyyy, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_zzzz, g_x_0_xxxxx_zzzzz, g_x_0_xxxxxz_xxxx, g_x_0_xxxxxz_xxxy, g_x_0_xxxxxz_xxxz, g_x_0_xxxxxz_xxyy, g_x_0_xxxxxz_xxyz, g_x_0_xxxxxz_xxzz, g_x_0_xxxxxz_xyyy, g_x_0_xxxxxz_xyyz, g_x_0_xxxxxz_xyzz, g_x_0_xxxxxz_xzzz, g_x_0_xxxxxz_yyyy, g_x_0_xxxxxz_yyyz, g_x_0_xxxxxz_yyzz, g_x_0_xxxxxz_yzzz, g_x_0_xxxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_xxxx[k] = -g_x_0_xxxxx_xxxx[k] * cd_z[k] + g_x_0_xxxxx_xxxxz[k];

                g_x_0_xxxxxz_xxxy[k] = -g_x_0_xxxxx_xxxy[k] * cd_z[k] + g_x_0_xxxxx_xxxyz[k];

                g_x_0_xxxxxz_xxxz[k] = -g_x_0_xxxxx_xxxz[k] * cd_z[k] + g_x_0_xxxxx_xxxzz[k];

                g_x_0_xxxxxz_xxyy[k] = -g_x_0_xxxxx_xxyy[k] * cd_z[k] + g_x_0_xxxxx_xxyyz[k];

                g_x_0_xxxxxz_xxyz[k] = -g_x_0_xxxxx_xxyz[k] * cd_z[k] + g_x_0_xxxxx_xxyzz[k];

                g_x_0_xxxxxz_xxzz[k] = -g_x_0_xxxxx_xxzz[k] * cd_z[k] + g_x_0_xxxxx_xxzzz[k];

                g_x_0_xxxxxz_xyyy[k] = -g_x_0_xxxxx_xyyy[k] * cd_z[k] + g_x_0_xxxxx_xyyyz[k];

                g_x_0_xxxxxz_xyyz[k] = -g_x_0_xxxxx_xyyz[k] * cd_z[k] + g_x_0_xxxxx_xyyzz[k];

                g_x_0_xxxxxz_xyzz[k] = -g_x_0_xxxxx_xyzz[k] * cd_z[k] + g_x_0_xxxxx_xyzzz[k];

                g_x_0_xxxxxz_xzzz[k] = -g_x_0_xxxxx_xzzz[k] * cd_z[k] + g_x_0_xxxxx_xzzzz[k];

                g_x_0_xxxxxz_yyyy[k] = -g_x_0_xxxxx_yyyy[k] * cd_z[k] + g_x_0_xxxxx_yyyyz[k];

                g_x_0_xxxxxz_yyyz[k] = -g_x_0_xxxxx_yyyz[k] * cd_z[k] + g_x_0_xxxxx_yyyzz[k];

                g_x_0_xxxxxz_yyzz[k] = -g_x_0_xxxxx_yyzz[k] * cd_z[k] + g_x_0_xxxxx_yyzzz[k];

                g_x_0_xxxxxz_yzzz[k] = -g_x_0_xxxxx_yzzz[k] * cd_z[k] + g_x_0_xxxxx_yzzzz[k];

                g_x_0_xxxxxz_zzzz[k] = -g_x_0_xxxxx_zzzz[k] * cd_z[k] + g_x_0_xxxxx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxxyy_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxxyy_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxxyy_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxxyy_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxxyy_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxxyy_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxxyy_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxxyy_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxxyy_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxxyy_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxxyy_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxxyy_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxxyy_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxxyy_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxy_xxxx, g_x_0_xxxxy_xxxxy, g_x_0_xxxxy_xxxy, g_x_0_xxxxy_xxxyy, g_x_0_xxxxy_xxxyz, g_x_0_xxxxy_xxxz, g_x_0_xxxxy_xxyy, g_x_0_xxxxy_xxyyy, g_x_0_xxxxy_xxyyz, g_x_0_xxxxy_xxyz, g_x_0_xxxxy_xxyzz, g_x_0_xxxxy_xxzz, g_x_0_xxxxy_xyyy, g_x_0_xxxxy_xyyyy, g_x_0_xxxxy_xyyyz, g_x_0_xxxxy_xyyz, g_x_0_xxxxy_xyyzz, g_x_0_xxxxy_xyzz, g_x_0_xxxxy_xyzzz, g_x_0_xxxxy_xzzz, g_x_0_xxxxy_yyyy, g_x_0_xxxxy_yyyyy, g_x_0_xxxxy_yyyyz, g_x_0_xxxxy_yyyz, g_x_0_xxxxy_yyyzz, g_x_0_xxxxy_yyzz, g_x_0_xxxxy_yyzzz, g_x_0_xxxxy_yzzz, g_x_0_xxxxy_yzzzz, g_x_0_xxxxy_zzzz, g_x_0_xxxxyy_xxxx, g_x_0_xxxxyy_xxxy, g_x_0_xxxxyy_xxxz, g_x_0_xxxxyy_xxyy, g_x_0_xxxxyy_xxyz, g_x_0_xxxxyy_xxzz, g_x_0_xxxxyy_xyyy, g_x_0_xxxxyy_xyyz, g_x_0_xxxxyy_xyzz, g_x_0_xxxxyy_xzzz, g_x_0_xxxxyy_yyyy, g_x_0_xxxxyy_yyyz, g_x_0_xxxxyy_yyzz, g_x_0_xxxxyy_yzzz, g_x_0_xxxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_xxxx[k] = -g_x_0_xxxxy_xxxx[k] * cd_y[k] + g_x_0_xxxxy_xxxxy[k];

                g_x_0_xxxxyy_xxxy[k] = -g_x_0_xxxxy_xxxy[k] * cd_y[k] + g_x_0_xxxxy_xxxyy[k];

                g_x_0_xxxxyy_xxxz[k] = -g_x_0_xxxxy_xxxz[k] * cd_y[k] + g_x_0_xxxxy_xxxyz[k];

                g_x_0_xxxxyy_xxyy[k] = -g_x_0_xxxxy_xxyy[k] * cd_y[k] + g_x_0_xxxxy_xxyyy[k];

                g_x_0_xxxxyy_xxyz[k] = -g_x_0_xxxxy_xxyz[k] * cd_y[k] + g_x_0_xxxxy_xxyyz[k];

                g_x_0_xxxxyy_xxzz[k] = -g_x_0_xxxxy_xxzz[k] * cd_y[k] + g_x_0_xxxxy_xxyzz[k];

                g_x_0_xxxxyy_xyyy[k] = -g_x_0_xxxxy_xyyy[k] * cd_y[k] + g_x_0_xxxxy_xyyyy[k];

                g_x_0_xxxxyy_xyyz[k] = -g_x_0_xxxxy_xyyz[k] * cd_y[k] + g_x_0_xxxxy_xyyyz[k];

                g_x_0_xxxxyy_xyzz[k] = -g_x_0_xxxxy_xyzz[k] * cd_y[k] + g_x_0_xxxxy_xyyzz[k];

                g_x_0_xxxxyy_xzzz[k] = -g_x_0_xxxxy_xzzz[k] * cd_y[k] + g_x_0_xxxxy_xyzzz[k];

                g_x_0_xxxxyy_yyyy[k] = -g_x_0_xxxxy_yyyy[k] * cd_y[k] + g_x_0_xxxxy_yyyyy[k];

                g_x_0_xxxxyy_yyyz[k] = -g_x_0_xxxxy_yyyz[k] * cd_y[k] + g_x_0_xxxxy_yyyyz[k];

                g_x_0_xxxxyy_yyzz[k] = -g_x_0_xxxxy_yyzz[k] * cd_y[k] + g_x_0_xxxxy_yyyzz[k];

                g_x_0_xxxxyy_yzzz[k] = -g_x_0_xxxxy_yzzz[k] * cd_y[k] + g_x_0_xxxxy_yyzzz[k];

                g_x_0_xxxxyy_zzzz[k] = -g_x_0_xxxxy_zzzz[k] * cd_y[k] + g_x_0_xxxxy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxxyz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxxyz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxxxyz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxxyz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxxyz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxxyz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxxyz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxxyz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxxyz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxxyz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxxyz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxxyz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxxyz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxxyz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 74);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxyz_xxxx, g_x_0_xxxxyz_xxxy, g_x_0_xxxxyz_xxxz, g_x_0_xxxxyz_xxyy, g_x_0_xxxxyz_xxyz, g_x_0_xxxxyz_xxzz, g_x_0_xxxxyz_xyyy, g_x_0_xxxxyz_xyyz, g_x_0_xxxxyz_xyzz, g_x_0_xxxxyz_xzzz, g_x_0_xxxxyz_yyyy, g_x_0_xxxxyz_yyyz, g_x_0_xxxxyz_yyzz, g_x_0_xxxxyz_yzzz, g_x_0_xxxxyz_zzzz, g_x_0_xxxxz_xxxx, g_x_0_xxxxz_xxxxy, g_x_0_xxxxz_xxxy, g_x_0_xxxxz_xxxyy, g_x_0_xxxxz_xxxyz, g_x_0_xxxxz_xxxz, g_x_0_xxxxz_xxyy, g_x_0_xxxxz_xxyyy, g_x_0_xxxxz_xxyyz, g_x_0_xxxxz_xxyz, g_x_0_xxxxz_xxyzz, g_x_0_xxxxz_xxzz, g_x_0_xxxxz_xyyy, g_x_0_xxxxz_xyyyy, g_x_0_xxxxz_xyyyz, g_x_0_xxxxz_xyyz, g_x_0_xxxxz_xyyzz, g_x_0_xxxxz_xyzz, g_x_0_xxxxz_xyzzz, g_x_0_xxxxz_xzzz, g_x_0_xxxxz_yyyy, g_x_0_xxxxz_yyyyy, g_x_0_xxxxz_yyyyz, g_x_0_xxxxz_yyyz, g_x_0_xxxxz_yyyzz, g_x_0_xxxxz_yyzz, g_x_0_xxxxz_yyzzz, g_x_0_xxxxz_yzzz, g_x_0_xxxxz_yzzzz, g_x_0_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_xxxx[k] = -g_x_0_xxxxz_xxxx[k] * cd_y[k] + g_x_0_xxxxz_xxxxy[k];

                g_x_0_xxxxyz_xxxy[k] = -g_x_0_xxxxz_xxxy[k] * cd_y[k] + g_x_0_xxxxz_xxxyy[k];

                g_x_0_xxxxyz_xxxz[k] = -g_x_0_xxxxz_xxxz[k] * cd_y[k] + g_x_0_xxxxz_xxxyz[k];

                g_x_0_xxxxyz_xxyy[k] = -g_x_0_xxxxz_xxyy[k] * cd_y[k] + g_x_0_xxxxz_xxyyy[k];

                g_x_0_xxxxyz_xxyz[k] = -g_x_0_xxxxz_xxyz[k] * cd_y[k] + g_x_0_xxxxz_xxyyz[k];

                g_x_0_xxxxyz_xxzz[k] = -g_x_0_xxxxz_xxzz[k] * cd_y[k] + g_x_0_xxxxz_xxyzz[k];

                g_x_0_xxxxyz_xyyy[k] = -g_x_0_xxxxz_xyyy[k] * cd_y[k] + g_x_0_xxxxz_xyyyy[k];

                g_x_0_xxxxyz_xyyz[k] = -g_x_0_xxxxz_xyyz[k] * cd_y[k] + g_x_0_xxxxz_xyyyz[k];

                g_x_0_xxxxyz_xyzz[k] = -g_x_0_xxxxz_xyzz[k] * cd_y[k] + g_x_0_xxxxz_xyyzz[k];

                g_x_0_xxxxyz_xzzz[k] = -g_x_0_xxxxz_xzzz[k] * cd_y[k] + g_x_0_xxxxz_xyzzz[k];

                g_x_0_xxxxyz_yyyy[k] = -g_x_0_xxxxz_yyyy[k] * cd_y[k] + g_x_0_xxxxz_yyyyy[k];

                g_x_0_xxxxyz_yyyz[k] = -g_x_0_xxxxz_yyyz[k] * cd_y[k] + g_x_0_xxxxz_yyyyz[k];

                g_x_0_xxxxyz_yyzz[k] = -g_x_0_xxxxz_yyzz[k] * cd_y[k] + g_x_0_xxxxz_yyyzz[k];

                g_x_0_xxxxyz_yzzz[k] = -g_x_0_xxxxz_yzzz[k] * cd_y[k] + g_x_0_xxxxz_yyzzz[k];

                g_x_0_xxxxyz_zzzz[k] = -g_x_0_xxxxz_zzzz[k] * cd_y[k] + g_x_0_xxxxz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxxzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxxzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxxzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxxzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxxzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxxzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxxzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxxzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxxxzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxxxzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxxxzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxxxzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxxxzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxxxzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxz_xxxx, g_x_0_xxxxz_xxxxz, g_x_0_xxxxz_xxxy, g_x_0_xxxxz_xxxyz, g_x_0_xxxxz_xxxz, g_x_0_xxxxz_xxxzz, g_x_0_xxxxz_xxyy, g_x_0_xxxxz_xxyyz, g_x_0_xxxxz_xxyz, g_x_0_xxxxz_xxyzz, g_x_0_xxxxz_xxzz, g_x_0_xxxxz_xxzzz, g_x_0_xxxxz_xyyy, g_x_0_xxxxz_xyyyz, g_x_0_xxxxz_xyyz, g_x_0_xxxxz_xyyzz, g_x_0_xxxxz_xyzz, g_x_0_xxxxz_xyzzz, g_x_0_xxxxz_xzzz, g_x_0_xxxxz_xzzzz, g_x_0_xxxxz_yyyy, g_x_0_xxxxz_yyyyz, g_x_0_xxxxz_yyyz, g_x_0_xxxxz_yyyzz, g_x_0_xxxxz_yyzz, g_x_0_xxxxz_yyzzz, g_x_0_xxxxz_yzzz, g_x_0_xxxxz_yzzzz, g_x_0_xxxxz_zzzz, g_x_0_xxxxz_zzzzz, g_x_0_xxxxzz_xxxx, g_x_0_xxxxzz_xxxy, g_x_0_xxxxzz_xxxz, g_x_0_xxxxzz_xxyy, g_x_0_xxxxzz_xxyz, g_x_0_xxxxzz_xxzz, g_x_0_xxxxzz_xyyy, g_x_0_xxxxzz_xyyz, g_x_0_xxxxzz_xyzz, g_x_0_xxxxzz_xzzz, g_x_0_xxxxzz_yyyy, g_x_0_xxxxzz_yyyz, g_x_0_xxxxzz_yyzz, g_x_0_xxxxzz_yzzz, g_x_0_xxxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_xxxx[k] = -g_x_0_xxxxz_xxxx[k] * cd_z[k] + g_x_0_xxxxz_xxxxz[k];

                g_x_0_xxxxzz_xxxy[k] = -g_x_0_xxxxz_xxxy[k] * cd_z[k] + g_x_0_xxxxz_xxxyz[k];

                g_x_0_xxxxzz_xxxz[k] = -g_x_0_xxxxz_xxxz[k] * cd_z[k] + g_x_0_xxxxz_xxxzz[k];

                g_x_0_xxxxzz_xxyy[k] = -g_x_0_xxxxz_xxyy[k] * cd_z[k] + g_x_0_xxxxz_xxyyz[k];

                g_x_0_xxxxzz_xxyz[k] = -g_x_0_xxxxz_xxyz[k] * cd_z[k] + g_x_0_xxxxz_xxyzz[k];

                g_x_0_xxxxzz_xxzz[k] = -g_x_0_xxxxz_xxzz[k] * cd_z[k] + g_x_0_xxxxz_xxzzz[k];

                g_x_0_xxxxzz_xyyy[k] = -g_x_0_xxxxz_xyyy[k] * cd_z[k] + g_x_0_xxxxz_xyyyz[k];

                g_x_0_xxxxzz_xyyz[k] = -g_x_0_xxxxz_xyyz[k] * cd_z[k] + g_x_0_xxxxz_xyyzz[k];

                g_x_0_xxxxzz_xyzz[k] = -g_x_0_xxxxz_xyzz[k] * cd_z[k] + g_x_0_xxxxz_xyzzz[k];

                g_x_0_xxxxzz_xzzz[k] = -g_x_0_xxxxz_xzzz[k] * cd_z[k] + g_x_0_xxxxz_xzzzz[k];

                g_x_0_xxxxzz_yyyy[k] = -g_x_0_xxxxz_yyyy[k] * cd_z[k] + g_x_0_xxxxz_yyyyz[k];

                g_x_0_xxxxzz_yyyz[k] = -g_x_0_xxxxz_yyyz[k] * cd_z[k] + g_x_0_xxxxz_yyyzz[k];

                g_x_0_xxxxzz_yyzz[k] = -g_x_0_xxxxz_yyzz[k] * cd_z[k] + g_x_0_xxxxz_yyzzz[k];

                g_x_0_xxxxzz_yzzz[k] = -g_x_0_xxxxz_yzzz[k] * cd_z[k] + g_x_0_xxxxz_yzzzz[k];

                g_x_0_xxxxzz_zzzz[k] = -g_x_0_xxxxz_zzzz[k] * cd_z[k] + g_x_0_xxxxz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxxyyy_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxxyyy_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxxyyy_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxxyyy_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxxyyy_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxxyyy_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxxyyy_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxxyyy_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxxyyy_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxxyyy_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxxyyy_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxxyyy_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxxyyy_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxxyyy_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyy_xxxx, g_x_0_xxxyy_xxxxy, g_x_0_xxxyy_xxxy, g_x_0_xxxyy_xxxyy, g_x_0_xxxyy_xxxyz, g_x_0_xxxyy_xxxz, g_x_0_xxxyy_xxyy, g_x_0_xxxyy_xxyyy, g_x_0_xxxyy_xxyyz, g_x_0_xxxyy_xxyz, g_x_0_xxxyy_xxyzz, g_x_0_xxxyy_xxzz, g_x_0_xxxyy_xyyy, g_x_0_xxxyy_xyyyy, g_x_0_xxxyy_xyyyz, g_x_0_xxxyy_xyyz, g_x_0_xxxyy_xyyzz, g_x_0_xxxyy_xyzz, g_x_0_xxxyy_xyzzz, g_x_0_xxxyy_xzzz, g_x_0_xxxyy_yyyy, g_x_0_xxxyy_yyyyy, g_x_0_xxxyy_yyyyz, g_x_0_xxxyy_yyyz, g_x_0_xxxyy_yyyzz, g_x_0_xxxyy_yyzz, g_x_0_xxxyy_yyzzz, g_x_0_xxxyy_yzzz, g_x_0_xxxyy_yzzzz, g_x_0_xxxyy_zzzz, g_x_0_xxxyyy_xxxx, g_x_0_xxxyyy_xxxy, g_x_0_xxxyyy_xxxz, g_x_0_xxxyyy_xxyy, g_x_0_xxxyyy_xxyz, g_x_0_xxxyyy_xxzz, g_x_0_xxxyyy_xyyy, g_x_0_xxxyyy_xyyz, g_x_0_xxxyyy_xyzz, g_x_0_xxxyyy_xzzz, g_x_0_xxxyyy_yyyy, g_x_0_xxxyyy_yyyz, g_x_0_xxxyyy_yyzz, g_x_0_xxxyyy_yzzz, g_x_0_xxxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_xxxx[k] = -g_x_0_xxxyy_xxxx[k] * cd_y[k] + g_x_0_xxxyy_xxxxy[k];

                g_x_0_xxxyyy_xxxy[k] = -g_x_0_xxxyy_xxxy[k] * cd_y[k] + g_x_0_xxxyy_xxxyy[k];

                g_x_0_xxxyyy_xxxz[k] = -g_x_0_xxxyy_xxxz[k] * cd_y[k] + g_x_0_xxxyy_xxxyz[k];

                g_x_0_xxxyyy_xxyy[k] = -g_x_0_xxxyy_xxyy[k] * cd_y[k] + g_x_0_xxxyy_xxyyy[k];

                g_x_0_xxxyyy_xxyz[k] = -g_x_0_xxxyy_xxyz[k] * cd_y[k] + g_x_0_xxxyy_xxyyz[k];

                g_x_0_xxxyyy_xxzz[k] = -g_x_0_xxxyy_xxzz[k] * cd_y[k] + g_x_0_xxxyy_xxyzz[k];

                g_x_0_xxxyyy_xyyy[k] = -g_x_0_xxxyy_xyyy[k] * cd_y[k] + g_x_0_xxxyy_xyyyy[k];

                g_x_0_xxxyyy_xyyz[k] = -g_x_0_xxxyy_xyyz[k] * cd_y[k] + g_x_0_xxxyy_xyyyz[k];

                g_x_0_xxxyyy_xyzz[k] = -g_x_0_xxxyy_xyzz[k] * cd_y[k] + g_x_0_xxxyy_xyyzz[k];

                g_x_0_xxxyyy_xzzz[k] = -g_x_0_xxxyy_xzzz[k] * cd_y[k] + g_x_0_xxxyy_xyzzz[k];

                g_x_0_xxxyyy_yyyy[k] = -g_x_0_xxxyy_yyyy[k] * cd_y[k] + g_x_0_xxxyy_yyyyy[k];

                g_x_0_xxxyyy_yyyz[k] = -g_x_0_xxxyy_yyyz[k] * cd_y[k] + g_x_0_xxxyy_yyyyz[k];

                g_x_0_xxxyyy_yyzz[k] = -g_x_0_xxxyy_yyzz[k] * cd_y[k] + g_x_0_xxxyy_yyyzz[k];

                g_x_0_xxxyyy_yzzz[k] = -g_x_0_xxxyy_yzzz[k] * cd_y[k] + g_x_0_xxxyy_yyzzz[k];

                g_x_0_xxxyyy_zzzz[k] = -g_x_0_xxxyy_zzzz[k] * cd_y[k] + g_x_0_xxxyy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxxyyz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxxyyz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxxyyz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxxyyz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxxyyz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxxyyz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxxyyz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxxyyz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxxyyz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxxyyz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxxyyz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxxyyz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxxyyz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxxyyz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyyz_xxxx, g_x_0_xxxyyz_xxxy, g_x_0_xxxyyz_xxxz, g_x_0_xxxyyz_xxyy, g_x_0_xxxyyz_xxyz, g_x_0_xxxyyz_xxzz, g_x_0_xxxyyz_xyyy, g_x_0_xxxyyz_xyyz, g_x_0_xxxyyz_xyzz, g_x_0_xxxyyz_xzzz, g_x_0_xxxyyz_yyyy, g_x_0_xxxyyz_yyyz, g_x_0_xxxyyz_yyzz, g_x_0_xxxyyz_yzzz, g_x_0_xxxyyz_zzzz, g_x_0_xxxyz_xxxx, g_x_0_xxxyz_xxxxy, g_x_0_xxxyz_xxxy, g_x_0_xxxyz_xxxyy, g_x_0_xxxyz_xxxyz, g_x_0_xxxyz_xxxz, g_x_0_xxxyz_xxyy, g_x_0_xxxyz_xxyyy, g_x_0_xxxyz_xxyyz, g_x_0_xxxyz_xxyz, g_x_0_xxxyz_xxyzz, g_x_0_xxxyz_xxzz, g_x_0_xxxyz_xyyy, g_x_0_xxxyz_xyyyy, g_x_0_xxxyz_xyyyz, g_x_0_xxxyz_xyyz, g_x_0_xxxyz_xyyzz, g_x_0_xxxyz_xyzz, g_x_0_xxxyz_xyzzz, g_x_0_xxxyz_xzzz, g_x_0_xxxyz_yyyy, g_x_0_xxxyz_yyyyy, g_x_0_xxxyz_yyyyz, g_x_0_xxxyz_yyyz, g_x_0_xxxyz_yyyzz, g_x_0_xxxyz_yyzz, g_x_0_xxxyz_yyzzz, g_x_0_xxxyz_yzzz, g_x_0_xxxyz_yzzzz, g_x_0_xxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_xxxx[k] = -g_x_0_xxxyz_xxxx[k] * cd_y[k] + g_x_0_xxxyz_xxxxy[k];

                g_x_0_xxxyyz_xxxy[k] = -g_x_0_xxxyz_xxxy[k] * cd_y[k] + g_x_0_xxxyz_xxxyy[k];

                g_x_0_xxxyyz_xxxz[k] = -g_x_0_xxxyz_xxxz[k] * cd_y[k] + g_x_0_xxxyz_xxxyz[k];

                g_x_0_xxxyyz_xxyy[k] = -g_x_0_xxxyz_xxyy[k] * cd_y[k] + g_x_0_xxxyz_xxyyy[k];

                g_x_0_xxxyyz_xxyz[k] = -g_x_0_xxxyz_xxyz[k] * cd_y[k] + g_x_0_xxxyz_xxyyz[k];

                g_x_0_xxxyyz_xxzz[k] = -g_x_0_xxxyz_xxzz[k] * cd_y[k] + g_x_0_xxxyz_xxyzz[k];

                g_x_0_xxxyyz_xyyy[k] = -g_x_0_xxxyz_xyyy[k] * cd_y[k] + g_x_0_xxxyz_xyyyy[k];

                g_x_0_xxxyyz_xyyz[k] = -g_x_0_xxxyz_xyyz[k] * cd_y[k] + g_x_0_xxxyz_xyyyz[k];

                g_x_0_xxxyyz_xyzz[k] = -g_x_0_xxxyz_xyzz[k] * cd_y[k] + g_x_0_xxxyz_xyyzz[k];

                g_x_0_xxxyyz_xzzz[k] = -g_x_0_xxxyz_xzzz[k] * cd_y[k] + g_x_0_xxxyz_xyzzz[k];

                g_x_0_xxxyyz_yyyy[k] = -g_x_0_xxxyz_yyyy[k] * cd_y[k] + g_x_0_xxxyz_yyyyy[k];

                g_x_0_xxxyyz_yyyz[k] = -g_x_0_xxxyz_yyyz[k] * cd_y[k] + g_x_0_xxxyz_yyyyz[k];

                g_x_0_xxxyyz_yyzz[k] = -g_x_0_xxxyz_yyzz[k] * cd_y[k] + g_x_0_xxxyz_yyyzz[k];

                g_x_0_xxxyyz_yzzz[k] = -g_x_0_xxxyz_yzzz[k] * cd_y[k] + g_x_0_xxxyz_yyzzz[k];

                g_x_0_xxxyyz_zzzz[k] = -g_x_0_xxxyz_zzzz[k] * cd_y[k] + g_x_0_xxxyz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxxyzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxxyzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxxyzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxxyzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxxyzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxxyzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxxyzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxxyzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxxyzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxxyzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxxyzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxxyzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxxyzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxxyzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 134);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyzz_xxxx, g_x_0_xxxyzz_xxxy, g_x_0_xxxyzz_xxxz, g_x_0_xxxyzz_xxyy, g_x_0_xxxyzz_xxyz, g_x_0_xxxyzz_xxzz, g_x_0_xxxyzz_xyyy, g_x_0_xxxyzz_xyyz, g_x_0_xxxyzz_xyzz, g_x_0_xxxyzz_xzzz, g_x_0_xxxyzz_yyyy, g_x_0_xxxyzz_yyyz, g_x_0_xxxyzz_yyzz, g_x_0_xxxyzz_yzzz, g_x_0_xxxyzz_zzzz, g_x_0_xxxzz_xxxx, g_x_0_xxxzz_xxxxy, g_x_0_xxxzz_xxxy, g_x_0_xxxzz_xxxyy, g_x_0_xxxzz_xxxyz, g_x_0_xxxzz_xxxz, g_x_0_xxxzz_xxyy, g_x_0_xxxzz_xxyyy, g_x_0_xxxzz_xxyyz, g_x_0_xxxzz_xxyz, g_x_0_xxxzz_xxyzz, g_x_0_xxxzz_xxzz, g_x_0_xxxzz_xyyy, g_x_0_xxxzz_xyyyy, g_x_0_xxxzz_xyyyz, g_x_0_xxxzz_xyyz, g_x_0_xxxzz_xyyzz, g_x_0_xxxzz_xyzz, g_x_0_xxxzz_xyzzz, g_x_0_xxxzz_xzzz, g_x_0_xxxzz_yyyy, g_x_0_xxxzz_yyyyy, g_x_0_xxxzz_yyyyz, g_x_0_xxxzz_yyyz, g_x_0_xxxzz_yyyzz, g_x_0_xxxzz_yyzz, g_x_0_xxxzz_yyzzz, g_x_0_xxxzz_yzzz, g_x_0_xxxzz_yzzzz, g_x_0_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_xxxx[k] = -g_x_0_xxxzz_xxxx[k] * cd_y[k] + g_x_0_xxxzz_xxxxy[k];

                g_x_0_xxxyzz_xxxy[k] = -g_x_0_xxxzz_xxxy[k] * cd_y[k] + g_x_0_xxxzz_xxxyy[k];

                g_x_0_xxxyzz_xxxz[k] = -g_x_0_xxxzz_xxxz[k] * cd_y[k] + g_x_0_xxxzz_xxxyz[k];

                g_x_0_xxxyzz_xxyy[k] = -g_x_0_xxxzz_xxyy[k] * cd_y[k] + g_x_0_xxxzz_xxyyy[k];

                g_x_0_xxxyzz_xxyz[k] = -g_x_0_xxxzz_xxyz[k] * cd_y[k] + g_x_0_xxxzz_xxyyz[k];

                g_x_0_xxxyzz_xxzz[k] = -g_x_0_xxxzz_xxzz[k] * cd_y[k] + g_x_0_xxxzz_xxyzz[k];

                g_x_0_xxxyzz_xyyy[k] = -g_x_0_xxxzz_xyyy[k] * cd_y[k] + g_x_0_xxxzz_xyyyy[k];

                g_x_0_xxxyzz_xyyz[k] = -g_x_0_xxxzz_xyyz[k] * cd_y[k] + g_x_0_xxxzz_xyyyz[k];

                g_x_0_xxxyzz_xyzz[k] = -g_x_0_xxxzz_xyzz[k] * cd_y[k] + g_x_0_xxxzz_xyyzz[k];

                g_x_0_xxxyzz_xzzz[k] = -g_x_0_xxxzz_xzzz[k] * cd_y[k] + g_x_0_xxxzz_xyzzz[k];

                g_x_0_xxxyzz_yyyy[k] = -g_x_0_xxxzz_yyyy[k] * cd_y[k] + g_x_0_xxxzz_yyyyy[k];

                g_x_0_xxxyzz_yyyz[k] = -g_x_0_xxxzz_yyyz[k] * cd_y[k] + g_x_0_xxxzz_yyyyz[k];

                g_x_0_xxxyzz_yyzz[k] = -g_x_0_xxxzz_yyzz[k] * cd_y[k] + g_x_0_xxxzz_yyyzz[k];

                g_x_0_xxxyzz_yzzz[k] = -g_x_0_xxxzz_yzzz[k] * cd_y[k] + g_x_0_xxxzz_yyzzz[k];

                g_x_0_xxxyzz_zzzz[k] = -g_x_0_xxxzz_zzzz[k] * cd_y[k] + g_x_0_xxxzz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxxzzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxxzzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxxzzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxxzzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xxxzzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxxzzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxxzzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxxzzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxxzzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxxzzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxxzzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxxzzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxxzzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxxzzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_z, g_x_0_xxxzz_xxxx, g_x_0_xxxzz_xxxxz, g_x_0_xxxzz_xxxy, g_x_0_xxxzz_xxxyz, g_x_0_xxxzz_xxxz, g_x_0_xxxzz_xxxzz, g_x_0_xxxzz_xxyy, g_x_0_xxxzz_xxyyz, g_x_0_xxxzz_xxyz, g_x_0_xxxzz_xxyzz, g_x_0_xxxzz_xxzz, g_x_0_xxxzz_xxzzz, g_x_0_xxxzz_xyyy, g_x_0_xxxzz_xyyyz, g_x_0_xxxzz_xyyz, g_x_0_xxxzz_xyyzz, g_x_0_xxxzz_xyzz, g_x_0_xxxzz_xyzzz, g_x_0_xxxzz_xzzz, g_x_0_xxxzz_xzzzz, g_x_0_xxxzz_yyyy, g_x_0_xxxzz_yyyyz, g_x_0_xxxzz_yyyz, g_x_0_xxxzz_yyyzz, g_x_0_xxxzz_yyzz, g_x_0_xxxzz_yyzzz, g_x_0_xxxzz_yzzz, g_x_0_xxxzz_yzzzz, g_x_0_xxxzz_zzzz, g_x_0_xxxzz_zzzzz, g_x_0_xxxzzz_xxxx, g_x_0_xxxzzz_xxxy, g_x_0_xxxzzz_xxxz, g_x_0_xxxzzz_xxyy, g_x_0_xxxzzz_xxyz, g_x_0_xxxzzz_xxzz, g_x_0_xxxzzz_xyyy, g_x_0_xxxzzz_xyyz, g_x_0_xxxzzz_xyzz, g_x_0_xxxzzz_xzzz, g_x_0_xxxzzz_yyyy, g_x_0_xxxzzz_yyyz, g_x_0_xxxzzz_yyzz, g_x_0_xxxzzz_yzzz, g_x_0_xxxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_xxxx[k] = -g_x_0_xxxzz_xxxx[k] * cd_z[k] + g_x_0_xxxzz_xxxxz[k];

                g_x_0_xxxzzz_xxxy[k] = -g_x_0_xxxzz_xxxy[k] * cd_z[k] + g_x_0_xxxzz_xxxyz[k];

                g_x_0_xxxzzz_xxxz[k] = -g_x_0_xxxzz_xxxz[k] * cd_z[k] + g_x_0_xxxzz_xxxzz[k];

                g_x_0_xxxzzz_xxyy[k] = -g_x_0_xxxzz_xxyy[k] * cd_z[k] + g_x_0_xxxzz_xxyyz[k];

                g_x_0_xxxzzz_xxyz[k] = -g_x_0_xxxzz_xxyz[k] * cd_z[k] + g_x_0_xxxzz_xxyzz[k];

                g_x_0_xxxzzz_xxzz[k] = -g_x_0_xxxzz_xxzz[k] * cd_z[k] + g_x_0_xxxzz_xxzzz[k];

                g_x_0_xxxzzz_xyyy[k] = -g_x_0_xxxzz_xyyy[k] * cd_z[k] + g_x_0_xxxzz_xyyyz[k];

                g_x_0_xxxzzz_xyyz[k] = -g_x_0_xxxzz_xyyz[k] * cd_z[k] + g_x_0_xxxzz_xyyzz[k];

                g_x_0_xxxzzz_xyzz[k] = -g_x_0_xxxzz_xyzz[k] * cd_z[k] + g_x_0_xxxzz_xyzzz[k];

                g_x_0_xxxzzz_xzzz[k] = -g_x_0_xxxzz_xzzz[k] * cd_z[k] + g_x_0_xxxzz_xzzzz[k];

                g_x_0_xxxzzz_yyyy[k] = -g_x_0_xxxzz_yyyy[k] * cd_z[k] + g_x_0_xxxzz_yyyyz[k];

                g_x_0_xxxzzz_yyyz[k] = -g_x_0_xxxzz_yyyz[k] * cd_z[k] + g_x_0_xxxzz_yyyzz[k];

                g_x_0_xxxzzz_yyzz[k] = -g_x_0_xxxzz_yyzz[k] * cd_z[k] + g_x_0_xxxzz_yyzzz[k];

                g_x_0_xxxzzz_yzzz[k] = -g_x_0_xxxzz_yzzz[k] * cd_z[k] + g_x_0_xxxzz_yzzzz[k];

                g_x_0_xxxzzz_zzzz[k] = -g_x_0_xxxzz_zzzz[k] * cd_z[k] + g_x_0_xxxzz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xxyyyy_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xxyyyy_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xxyyyy_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xxyyyy_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xxyyyy_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xxyyyy_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xxyyyy_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xxyyyy_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xxyyyy_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xxyyyy_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xxyyyy_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xxyyyy_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xxyyyy_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xxyyyy_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 164);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyy_xxxx, g_x_0_xxyyy_xxxxy, g_x_0_xxyyy_xxxy, g_x_0_xxyyy_xxxyy, g_x_0_xxyyy_xxxyz, g_x_0_xxyyy_xxxz, g_x_0_xxyyy_xxyy, g_x_0_xxyyy_xxyyy, g_x_0_xxyyy_xxyyz, g_x_0_xxyyy_xxyz, g_x_0_xxyyy_xxyzz, g_x_0_xxyyy_xxzz, g_x_0_xxyyy_xyyy, g_x_0_xxyyy_xyyyy, g_x_0_xxyyy_xyyyz, g_x_0_xxyyy_xyyz, g_x_0_xxyyy_xyyzz, g_x_0_xxyyy_xyzz, g_x_0_xxyyy_xyzzz, g_x_0_xxyyy_xzzz, g_x_0_xxyyy_yyyy, g_x_0_xxyyy_yyyyy, g_x_0_xxyyy_yyyyz, g_x_0_xxyyy_yyyz, g_x_0_xxyyy_yyyzz, g_x_0_xxyyy_yyzz, g_x_0_xxyyy_yyzzz, g_x_0_xxyyy_yzzz, g_x_0_xxyyy_yzzzz, g_x_0_xxyyy_zzzz, g_x_0_xxyyyy_xxxx, g_x_0_xxyyyy_xxxy, g_x_0_xxyyyy_xxxz, g_x_0_xxyyyy_xxyy, g_x_0_xxyyyy_xxyz, g_x_0_xxyyyy_xxzz, g_x_0_xxyyyy_xyyy, g_x_0_xxyyyy_xyyz, g_x_0_xxyyyy_xyzz, g_x_0_xxyyyy_xzzz, g_x_0_xxyyyy_yyyy, g_x_0_xxyyyy_yyyz, g_x_0_xxyyyy_yyzz, g_x_0_xxyyyy_yzzz, g_x_0_xxyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_xxxx[k] = -g_x_0_xxyyy_xxxx[k] * cd_y[k] + g_x_0_xxyyy_xxxxy[k];

                g_x_0_xxyyyy_xxxy[k] = -g_x_0_xxyyy_xxxy[k] * cd_y[k] + g_x_0_xxyyy_xxxyy[k];

                g_x_0_xxyyyy_xxxz[k] = -g_x_0_xxyyy_xxxz[k] * cd_y[k] + g_x_0_xxyyy_xxxyz[k];

                g_x_0_xxyyyy_xxyy[k] = -g_x_0_xxyyy_xxyy[k] * cd_y[k] + g_x_0_xxyyy_xxyyy[k];

                g_x_0_xxyyyy_xxyz[k] = -g_x_0_xxyyy_xxyz[k] * cd_y[k] + g_x_0_xxyyy_xxyyz[k];

                g_x_0_xxyyyy_xxzz[k] = -g_x_0_xxyyy_xxzz[k] * cd_y[k] + g_x_0_xxyyy_xxyzz[k];

                g_x_0_xxyyyy_xyyy[k] = -g_x_0_xxyyy_xyyy[k] * cd_y[k] + g_x_0_xxyyy_xyyyy[k];

                g_x_0_xxyyyy_xyyz[k] = -g_x_0_xxyyy_xyyz[k] * cd_y[k] + g_x_0_xxyyy_xyyyz[k];

                g_x_0_xxyyyy_xyzz[k] = -g_x_0_xxyyy_xyzz[k] * cd_y[k] + g_x_0_xxyyy_xyyzz[k];

                g_x_0_xxyyyy_xzzz[k] = -g_x_0_xxyyy_xzzz[k] * cd_y[k] + g_x_0_xxyyy_xyzzz[k];

                g_x_0_xxyyyy_yyyy[k] = -g_x_0_xxyyy_yyyy[k] * cd_y[k] + g_x_0_xxyyy_yyyyy[k];

                g_x_0_xxyyyy_yyyz[k] = -g_x_0_xxyyy_yyyz[k] * cd_y[k] + g_x_0_xxyyy_yyyyz[k];

                g_x_0_xxyyyy_yyzz[k] = -g_x_0_xxyyy_yyzz[k] * cd_y[k] + g_x_0_xxyyy_yyyzz[k];

                g_x_0_xxyyyy_yzzz[k] = -g_x_0_xxyyy_yzzz[k] * cd_y[k] + g_x_0_xxyyy_yyzzz[k];

                g_x_0_xxyyyy_zzzz[k] = -g_x_0_xxyyy_zzzz[k] * cd_y[k] + g_x_0_xxyyy_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xxyyyz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xxyyyz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xxyyyz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xxyyyz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xxyyyz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xxyyyz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xxyyyz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xxyyyz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xxyyyz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xxyyyz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xxyyyz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xxyyyz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xxyyyz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xxyyyz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyyz_xxxx, g_x_0_xxyyyz_xxxy, g_x_0_xxyyyz_xxxz, g_x_0_xxyyyz_xxyy, g_x_0_xxyyyz_xxyz, g_x_0_xxyyyz_xxzz, g_x_0_xxyyyz_xyyy, g_x_0_xxyyyz_xyyz, g_x_0_xxyyyz_xyzz, g_x_0_xxyyyz_xzzz, g_x_0_xxyyyz_yyyy, g_x_0_xxyyyz_yyyz, g_x_0_xxyyyz_yyzz, g_x_0_xxyyyz_yzzz, g_x_0_xxyyyz_zzzz, g_x_0_xxyyz_xxxx, g_x_0_xxyyz_xxxxy, g_x_0_xxyyz_xxxy, g_x_0_xxyyz_xxxyy, g_x_0_xxyyz_xxxyz, g_x_0_xxyyz_xxxz, g_x_0_xxyyz_xxyy, g_x_0_xxyyz_xxyyy, g_x_0_xxyyz_xxyyz, g_x_0_xxyyz_xxyz, g_x_0_xxyyz_xxyzz, g_x_0_xxyyz_xxzz, g_x_0_xxyyz_xyyy, g_x_0_xxyyz_xyyyy, g_x_0_xxyyz_xyyyz, g_x_0_xxyyz_xyyz, g_x_0_xxyyz_xyyzz, g_x_0_xxyyz_xyzz, g_x_0_xxyyz_xyzzz, g_x_0_xxyyz_xzzz, g_x_0_xxyyz_yyyy, g_x_0_xxyyz_yyyyy, g_x_0_xxyyz_yyyyz, g_x_0_xxyyz_yyyz, g_x_0_xxyyz_yyyzz, g_x_0_xxyyz_yyzz, g_x_0_xxyyz_yyzzz, g_x_0_xxyyz_yzzz, g_x_0_xxyyz_yzzzz, g_x_0_xxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_xxxx[k] = -g_x_0_xxyyz_xxxx[k] * cd_y[k] + g_x_0_xxyyz_xxxxy[k];

                g_x_0_xxyyyz_xxxy[k] = -g_x_0_xxyyz_xxxy[k] * cd_y[k] + g_x_0_xxyyz_xxxyy[k];

                g_x_0_xxyyyz_xxxz[k] = -g_x_0_xxyyz_xxxz[k] * cd_y[k] + g_x_0_xxyyz_xxxyz[k];

                g_x_0_xxyyyz_xxyy[k] = -g_x_0_xxyyz_xxyy[k] * cd_y[k] + g_x_0_xxyyz_xxyyy[k];

                g_x_0_xxyyyz_xxyz[k] = -g_x_0_xxyyz_xxyz[k] * cd_y[k] + g_x_0_xxyyz_xxyyz[k];

                g_x_0_xxyyyz_xxzz[k] = -g_x_0_xxyyz_xxzz[k] * cd_y[k] + g_x_0_xxyyz_xxyzz[k];

                g_x_0_xxyyyz_xyyy[k] = -g_x_0_xxyyz_xyyy[k] * cd_y[k] + g_x_0_xxyyz_xyyyy[k];

                g_x_0_xxyyyz_xyyz[k] = -g_x_0_xxyyz_xyyz[k] * cd_y[k] + g_x_0_xxyyz_xyyyz[k];

                g_x_0_xxyyyz_xyzz[k] = -g_x_0_xxyyz_xyzz[k] * cd_y[k] + g_x_0_xxyyz_xyyzz[k];

                g_x_0_xxyyyz_xzzz[k] = -g_x_0_xxyyz_xzzz[k] * cd_y[k] + g_x_0_xxyyz_xyzzz[k];

                g_x_0_xxyyyz_yyyy[k] = -g_x_0_xxyyz_yyyy[k] * cd_y[k] + g_x_0_xxyyz_yyyyy[k];

                g_x_0_xxyyyz_yyyz[k] = -g_x_0_xxyyz_yyyz[k] * cd_y[k] + g_x_0_xxyyz_yyyyz[k];

                g_x_0_xxyyyz_yyzz[k] = -g_x_0_xxyyz_yyzz[k] * cd_y[k] + g_x_0_xxyyz_yyyzz[k];

                g_x_0_xxyyyz_yzzz[k] = -g_x_0_xxyyz_yzzz[k] * cd_y[k] + g_x_0_xxyyz_yyzzz[k];

                g_x_0_xxyyyz_zzzz[k] = -g_x_0_xxyyz_zzzz[k] * cd_y[k] + g_x_0_xxyyz_yzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xxyyzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xxyyzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xxyyzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xxyyzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xxyyzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xxyyzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xxyyzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xxyyzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xxyyzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xxyyzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xxyyzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xxyyzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xxyyzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xxyyzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 194);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyzz_xxxx, g_x_0_xxyyzz_xxxy, g_x_0_xxyyzz_xxxz, g_x_0_xxyyzz_xxyy, g_x_0_xxyyzz_xxyz, g_x_0_xxyyzz_xxzz, g_x_0_xxyyzz_xyyy, g_x_0_xxyyzz_xyyz, g_x_0_xxyyzz_xyzz, g_x_0_xxyyzz_xzzz, g_x_0_xxyyzz_yyyy, g_x_0_xxyyzz_yyyz, g_x_0_xxyyzz_yyzz, g_x_0_xxyyzz_yzzz, g_x_0_xxyyzz_zzzz, g_x_0_xxyzz_xxxx, g_x_0_xxyzz_xxxxy, g_x_0_xxyzz_xxxy, g_x_0_xxyzz_xxxyy, g_x_0_xxyzz_xxxyz, g_x_0_xxyzz_xxxz, g_x_0_xxyzz_xxyy, g_x_0_xxyzz_xxyyy, g_x_0_xxyzz_xxyyz, g_x_0_xxyzz_xxyz, g_x_0_xxyzz_xxyzz, g_x_0_xxyzz_xxzz, g_x_0_xxyzz_xyyy, g_x_0_xxyzz_xyyyy, g_x_0_xxyzz_xyyyz, g_x_0_xxyzz_xyyz, g_x_0_xxyzz_xyyzz, g_x_0_xxyzz_xyzz, g_x_0_xxyzz_xyzzz, g_x_0_xxyzz_xzzz, g_x_0_xxyzz_yyyy, g_x_0_xxyzz_yyyyy, g_x_0_xxyzz_yyyyz, g_x_0_xxyzz_yyyz, g_x_0_xxyzz_yyyzz, g_x_0_xxyzz_yyzz, g_x_0_xxyzz_yyzzz, g_x_0_xxyzz_yzzz, g_x_0_xxyzz_yzzzz, g_x_0_xxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_xxxx[k] = -g_x_0_xxyzz_xxxx[k] * cd_y[k] + g_x_0_xxyzz_xxxxy[k];

                g_x_0_xxyyzz_xxxy[k] = -g_x_0_xxyzz_xxxy[k] * cd_y[k] + g_x_0_xxyzz_xxxyy[k];

                g_x_0_xxyyzz_xxxz[k] = -g_x_0_xxyzz_xxxz[k] * cd_y[k] + g_x_0_xxyzz_xxxyz[k];

                g_x_0_xxyyzz_xxyy[k] = -g_x_0_xxyzz_xxyy[k] * cd_y[k] + g_x_0_xxyzz_xxyyy[k];

                g_x_0_xxyyzz_xxyz[k] = -g_x_0_xxyzz_xxyz[k] * cd_y[k] + g_x_0_xxyzz_xxyyz[k];

                g_x_0_xxyyzz_xxzz[k] = -g_x_0_xxyzz_xxzz[k] * cd_y[k] + g_x_0_xxyzz_xxyzz[k];

                g_x_0_xxyyzz_xyyy[k] = -g_x_0_xxyzz_xyyy[k] * cd_y[k] + g_x_0_xxyzz_xyyyy[k];

                g_x_0_xxyyzz_xyyz[k] = -g_x_0_xxyzz_xyyz[k] * cd_y[k] + g_x_0_xxyzz_xyyyz[k];

                g_x_0_xxyyzz_xyzz[k] = -g_x_0_xxyzz_xyzz[k] * cd_y[k] + g_x_0_xxyzz_xyyzz[k];

                g_x_0_xxyyzz_xzzz[k] = -g_x_0_xxyzz_xzzz[k] * cd_y[k] + g_x_0_xxyzz_xyzzz[k];

                g_x_0_xxyyzz_yyyy[k] = -g_x_0_xxyzz_yyyy[k] * cd_y[k] + g_x_0_xxyzz_yyyyy[k];

                g_x_0_xxyyzz_yyyz[k] = -g_x_0_xxyzz_yyyz[k] * cd_y[k] + g_x_0_xxyzz_yyyyz[k];

                g_x_0_xxyyzz_yyzz[k] = -g_x_0_xxyzz_yyzz[k] * cd_y[k] + g_x_0_xxyzz_yyyzz[k];

                g_x_0_xxyyzz_yzzz[k] = -g_x_0_xxyzz_yzzz[k] * cd_y[k] + g_x_0_xxyzz_yyzzz[k];

                g_x_0_xxyyzz_zzzz[k] = -g_x_0_xxyzz_zzzz[k] * cd_y[k] + g_x_0_xxyzz_yzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xxyzzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xxyzzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xxyzzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xxyzzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xxyzzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xxyzzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xxyzzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xxyzzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xxyzzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xxyzzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xxyzzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xxyzzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xxyzzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xxyzzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzzz_xxxx, g_x_0_xxyzzz_xxxy, g_x_0_xxyzzz_xxxz, g_x_0_xxyzzz_xxyy, g_x_0_xxyzzz_xxyz, g_x_0_xxyzzz_xxzz, g_x_0_xxyzzz_xyyy, g_x_0_xxyzzz_xyyz, g_x_0_xxyzzz_xyzz, g_x_0_xxyzzz_xzzz, g_x_0_xxyzzz_yyyy, g_x_0_xxyzzz_yyyz, g_x_0_xxyzzz_yyzz, g_x_0_xxyzzz_yzzz, g_x_0_xxyzzz_zzzz, g_x_0_xxzzz_xxxx, g_x_0_xxzzz_xxxxy, g_x_0_xxzzz_xxxy, g_x_0_xxzzz_xxxyy, g_x_0_xxzzz_xxxyz, g_x_0_xxzzz_xxxz, g_x_0_xxzzz_xxyy, g_x_0_xxzzz_xxyyy, g_x_0_xxzzz_xxyyz, g_x_0_xxzzz_xxyz, g_x_0_xxzzz_xxyzz, g_x_0_xxzzz_xxzz, g_x_0_xxzzz_xyyy, g_x_0_xxzzz_xyyyy, g_x_0_xxzzz_xyyyz, g_x_0_xxzzz_xyyz, g_x_0_xxzzz_xyyzz, g_x_0_xxzzz_xyzz, g_x_0_xxzzz_xyzzz, g_x_0_xxzzz_xzzz, g_x_0_xxzzz_yyyy, g_x_0_xxzzz_yyyyy, g_x_0_xxzzz_yyyyz, g_x_0_xxzzz_yyyz, g_x_0_xxzzz_yyyzz, g_x_0_xxzzz_yyzz, g_x_0_xxzzz_yyzzz, g_x_0_xxzzz_yzzz, g_x_0_xxzzz_yzzzz, g_x_0_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_xxxx[k] = -g_x_0_xxzzz_xxxx[k] * cd_y[k] + g_x_0_xxzzz_xxxxy[k];

                g_x_0_xxyzzz_xxxy[k] = -g_x_0_xxzzz_xxxy[k] * cd_y[k] + g_x_0_xxzzz_xxxyy[k];

                g_x_0_xxyzzz_xxxz[k] = -g_x_0_xxzzz_xxxz[k] * cd_y[k] + g_x_0_xxzzz_xxxyz[k];

                g_x_0_xxyzzz_xxyy[k] = -g_x_0_xxzzz_xxyy[k] * cd_y[k] + g_x_0_xxzzz_xxyyy[k];

                g_x_0_xxyzzz_xxyz[k] = -g_x_0_xxzzz_xxyz[k] * cd_y[k] + g_x_0_xxzzz_xxyyz[k];

                g_x_0_xxyzzz_xxzz[k] = -g_x_0_xxzzz_xxzz[k] * cd_y[k] + g_x_0_xxzzz_xxyzz[k];

                g_x_0_xxyzzz_xyyy[k] = -g_x_0_xxzzz_xyyy[k] * cd_y[k] + g_x_0_xxzzz_xyyyy[k];

                g_x_0_xxyzzz_xyyz[k] = -g_x_0_xxzzz_xyyz[k] * cd_y[k] + g_x_0_xxzzz_xyyyz[k];

                g_x_0_xxyzzz_xyzz[k] = -g_x_0_xxzzz_xyzz[k] * cd_y[k] + g_x_0_xxzzz_xyyzz[k];

                g_x_0_xxyzzz_xzzz[k] = -g_x_0_xxzzz_xzzz[k] * cd_y[k] + g_x_0_xxzzz_xyzzz[k];

                g_x_0_xxyzzz_yyyy[k] = -g_x_0_xxzzz_yyyy[k] * cd_y[k] + g_x_0_xxzzz_yyyyy[k];

                g_x_0_xxyzzz_yyyz[k] = -g_x_0_xxzzz_yyyz[k] * cd_y[k] + g_x_0_xxzzz_yyyyz[k];

                g_x_0_xxyzzz_yyzz[k] = -g_x_0_xxzzz_yyzz[k] * cd_y[k] + g_x_0_xxzzz_yyyzz[k];

                g_x_0_xxyzzz_yzzz[k] = -g_x_0_xxzzz_yzzz[k] * cd_y[k] + g_x_0_xxzzz_yyzzz[k];

                g_x_0_xxyzzz_zzzz[k] = -g_x_0_xxzzz_zzzz[k] * cd_y[k] + g_x_0_xxzzz_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xxzzzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xxzzzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xxzzzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xxzzzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xxzzzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xxzzzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xxzzzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xxzzzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xxzzzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xxzzzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xxzzzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xxzzzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xxzzzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_xxzzzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 224);

            #pragma omp simd aligned(cd_z, g_x_0_xxzzz_xxxx, g_x_0_xxzzz_xxxxz, g_x_0_xxzzz_xxxy, g_x_0_xxzzz_xxxyz, g_x_0_xxzzz_xxxz, g_x_0_xxzzz_xxxzz, g_x_0_xxzzz_xxyy, g_x_0_xxzzz_xxyyz, g_x_0_xxzzz_xxyz, g_x_0_xxzzz_xxyzz, g_x_0_xxzzz_xxzz, g_x_0_xxzzz_xxzzz, g_x_0_xxzzz_xyyy, g_x_0_xxzzz_xyyyz, g_x_0_xxzzz_xyyz, g_x_0_xxzzz_xyyzz, g_x_0_xxzzz_xyzz, g_x_0_xxzzz_xyzzz, g_x_0_xxzzz_xzzz, g_x_0_xxzzz_xzzzz, g_x_0_xxzzz_yyyy, g_x_0_xxzzz_yyyyz, g_x_0_xxzzz_yyyz, g_x_0_xxzzz_yyyzz, g_x_0_xxzzz_yyzz, g_x_0_xxzzz_yyzzz, g_x_0_xxzzz_yzzz, g_x_0_xxzzz_yzzzz, g_x_0_xxzzz_zzzz, g_x_0_xxzzz_zzzzz, g_x_0_xxzzzz_xxxx, g_x_0_xxzzzz_xxxy, g_x_0_xxzzzz_xxxz, g_x_0_xxzzzz_xxyy, g_x_0_xxzzzz_xxyz, g_x_0_xxzzzz_xxzz, g_x_0_xxzzzz_xyyy, g_x_0_xxzzzz_xyyz, g_x_0_xxzzzz_xyzz, g_x_0_xxzzzz_xzzz, g_x_0_xxzzzz_yyyy, g_x_0_xxzzzz_yyyz, g_x_0_xxzzzz_yyzz, g_x_0_xxzzzz_yzzz, g_x_0_xxzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_xxxx[k] = -g_x_0_xxzzz_xxxx[k] * cd_z[k] + g_x_0_xxzzz_xxxxz[k];

                g_x_0_xxzzzz_xxxy[k] = -g_x_0_xxzzz_xxxy[k] * cd_z[k] + g_x_0_xxzzz_xxxyz[k];

                g_x_0_xxzzzz_xxxz[k] = -g_x_0_xxzzz_xxxz[k] * cd_z[k] + g_x_0_xxzzz_xxxzz[k];

                g_x_0_xxzzzz_xxyy[k] = -g_x_0_xxzzz_xxyy[k] * cd_z[k] + g_x_0_xxzzz_xxyyz[k];

                g_x_0_xxzzzz_xxyz[k] = -g_x_0_xxzzz_xxyz[k] * cd_z[k] + g_x_0_xxzzz_xxyzz[k];

                g_x_0_xxzzzz_xxzz[k] = -g_x_0_xxzzz_xxzz[k] * cd_z[k] + g_x_0_xxzzz_xxzzz[k];

                g_x_0_xxzzzz_xyyy[k] = -g_x_0_xxzzz_xyyy[k] * cd_z[k] + g_x_0_xxzzz_xyyyz[k];

                g_x_0_xxzzzz_xyyz[k] = -g_x_0_xxzzz_xyyz[k] * cd_z[k] + g_x_0_xxzzz_xyyzz[k];

                g_x_0_xxzzzz_xyzz[k] = -g_x_0_xxzzz_xyzz[k] * cd_z[k] + g_x_0_xxzzz_xyzzz[k];

                g_x_0_xxzzzz_xzzz[k] = -g_x_0_xxzzz_xzzz[k] * cd_z[k] + g_x_0_xxzzz_xzzzz[k];

                g_x_0_xxzzzz_yyyy[k] = -g_x_0_xxzzz_yyyy[k] * cd_z[k] + g_x_0_xxzzz_yyyyz[k];

                g_x_0_xxzzzz_yyyz[k] = -g_x_0_xxzzz_yyyz[k] * cd_z[k] + g_x_0_xxzzz_yyyzz[k];

                g_x_0_xxzzzz_yyzz[k] = -g_x_0_xxzzz_yyzz[k] * cd_z[k] + g_x_0_xxzzz_yyzzz[k];

                g_x_0_xxzzzz_yzzz[k] = -g_x_0_xxzzz_yzzz[k] * cd_z[k] + g_x_0_xxzzz_yzzzz[k];

                g_x_0_xxzzzz_zzzz[k] = -g_x_0_xxzzz_zzzz[k] * cd_z[k] + g_x_0_xxzzz_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_xyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_xyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_xyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_xyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_xyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_xyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_xyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_xyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_xyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_xyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_xyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_xyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_xyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_xyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 239);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyy_xxxx, g_x_0_xyyyy_xxxxy, g_x_0_xyyyy_xxxy, g_x_0_xyyyy_xxxyy, g_x_0_xyyyy_xxxyz, g_x_0_xyyyy_xxxz, g_x_0_xyyyy_xxyy, g_x_0_xyyyy_xxyyy, g_x_0_xyyyy_xxyyz, g_x_0_xyyyy_xxyz, g_x_0_xyyyy_xxyzz, g_x_0_xyyyy_xxzz, g_x_0_xyyyy_xyyy, g_x_0_xyyyy_xyyyy, g_x_0_xyyyy_xyyyz, g_x_0_xyyyy_xyyz, g_x_0_xyyyy_xyyzz, g_x_0_xyyyy_xyzz, g_x_0_xyyyy_xyzzz, g_x_0_xyyyy_xzzz, g_x_0_xyyyy_yyyy, g_x_0_xyyyy_yyyyy, g_x_0_xyyyy_yyyyz, g_x_0_xyyyy_yyyz, g_x_0_xyyyy_yyyzz, g_x_0_xyyyy_yyzz, g_x_0_xyyyy_yyzzz, g_x_0_xyyyy_yzzz, g_x_0_xyyyy_yzzzz, g_x_0_xyyyy_zzzz, g_x_0_xyyyyy_xxxx, g_x_0_xyyyyy_xxxy, g_x_0_xyyyyy_xxxz, g_x_0_xyyyyy_xxyy, g_x_0_xyyyyy_xxyz, g_x_0_xyyyyy_xxzz, g_x_0_xyyyyy_xyyy, g_x_0_xyyyyy_xyyz, g_x_0_xyyyyy_xyzz, g_x_0_xyyyyy_xzzz, g_x_0_xyyyyy_yyyy, g_x_0_xyyyyy_yyyz, g_x_0_xyyyyy_yyzz, g_x_0_xyyyyy_yzzz, g_x_0_xyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_xxxx[k] = -g_x_0_xyyyy_xxxx[k] * cd_y[k] + g_x_0_xyyyy_xxxxy[k];

                g_x_0_xyyyyy_xxxy[k] = -g_x_0_xyyyy_xxxy[k] * cd_y[k] + g_x_0_xyyyy_xxxyy[k];

                g_x_0_xyyyyy_xxxz[k] = -g_x_0_xyyyy_xxxz[k] * cd_y[k] + g_x_0_xyyyy_xxxyz[k];

                g_x_0_xyyyyy_xxyy[k] = -g_x_0_xyyyy_xxyy[k] * cd_y[k] + g_x_0_xyyyy_xxyyy[k];

                g_x_0_xyyyyy_xxyz[k] = -g_x_0_xyyyy_xxyz[k] * cd_y[k] + g_x_0_xyyyy_xxyyz[k];

                g_x_0_xyyyyy_xxzz[k] = -g_x_0_xyyyy_xxzz[k] * cd_y[k] + g_x_0_xyyyy_xxyzz[k];

                g_x_0_xyyyyy_xyyy[k] = -g_x_0_xyyyy_xyyy[k] * cd_y[k] + g_x_0_xyyyy_xyyyy[k];

                g_x_0_xyyyyy_xyyz[k] = -g_x_0_xyyyy_xyyz[k] * cd_y[k] + g_x_0_xyyyy_xyyyz[k];

                g_x_0_xyyyyy_xyzz[k] = -g_x_0_xyyyy_xyzz[k] * cd_y[k] + g_x_0_xyyyy_xyyzz[k];

                g_x_0_xyyyyy_xzzz[k] = -g_x_0_xyyyy_xzzz[k] * cd_y[k] + g_x_0_xyyyy_xyzzz[k];

                g_x_0_xyyyyy_yyyy[k] = -g_x_0_xyyyy_yyyy[k] * cd_y[k] + g_x_0_xyyyy_yyyyy[k];

                g_x_0_xyyyyy_yyyz[k] = -g_x_0_xyyyy_yyyz[k] * cd_y[k] + g_x_0_xyyyy_yyyyz[k];

                g_x_0_xyyyyy_yyzz[k] = -g_x_0_xyyyy_yyzz[k] * cd_y[k] + g_x_0_xyyyy_yyyzz[k];

                g_x_0_xyyyyy_yzzz[k] = -g_x_0_xyyyy_yzzz[k] * cd_y[k] + g_x_0_xyyyy_yyzzz[k];

                g_x_0_xyyyyy_zzzz[k] = -g_x_0_xyyyy_zzzz[k] * cd_y[k] + g_x_0_xyyyy_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_xyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_xyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_xyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_xyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_xyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_xyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_xyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_xyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_xyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_xyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_xyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_xyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_xyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_xyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 254);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyyz_xxxx, g_x_0_xyyyyz_xxxy, g_x_0_xyyyyz_xxxz, g_x_0_xyyyyz_xxyy, g_x_0_xyyyyz_xxyz, g_x_0_xyyyyz_xxzz, g_x_0_xyyyyz_xyyy, g_x_0_xyyyyz_xyyz, g_x_0_xyyyyz_xyzz, g_x_0_xyyyyz_xzzz, g_x_0_xyyyyz_yyyy, g_x_0_xyyyyz_yyyz, g_x_0_xyyyyz_yyzz, g_x_0_xyyyyz_yzzz, g_x_0_xyyyyz_zzzz, g_x_0_xyyyz_xxxx, g_x_0_xyyyz_xxxxy, g_x_0_xyyyz_xxxy, g_x_0_xyyyz_xxxyy, g_x_0_xyyyz_xxxyz, g_x_0_xyyyz_xxxz, g_x_0_xyyyz_xxyy, g_x_0_xyyyz_xxyyy, g_x_0_xyyyz_xxyyz, g_x_0_xyyyz_xxyz, g_x_0_xyyyz_xxyzz, g_x_0_xyyyz_xxzz, g_x_0_xyyyz_xyyy, g_x_0_xyyyz_xyyyy, g_x_0_xyyyz_xyyyz, g_x_0_xyyyz_xyyz, g_x_0_xyyyz_xyyzz, g_x_0_xyyyz_xyzz, g_x_0_xyyyz_xyzzz, g_x_0_xyyyz_xzzz, g_x_0_xyyyz_yyyy, g_x_0_xyyyz_yyyyy, g_x_0_xyyyz_yyyyz, g_x_0_xyyyz_yyyz, g_x_0_xyyyz_yyyzz, g_x_0_xyyyz_yyzz, g_x_0_xyyyz_yyzzz, g_x_0_xyyyz_yzzz, g_x_0_xyyyz_yzzzz, g_x_0_xyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_xxxx[k] = -g_x_0_xyyyz_xxxx[k] * cd_y[k] + g_x_0_xyyyz_xxxxy[k];

                g_x_0_xyyyyz_xxxy[k] = -g_x_0_xyyyz_xxxy[k] * cd_y[k] + g_x_0_xyyyz_xxxyy[k];

                g_x_0_xyyyyz_xxxz[k] = -g_x_0_xyyyz_xxxz[k] * cd_y[k] + g_x_0_xyyyz_xxxyz[k];

                g_x_0_xyyyyz_xxyy[k] = -g_x_0_xyyyz_xxyy[k] * cd_y[k] + g_x_0_xyyyz_xxyyy[k];

                g_x_0_xyyyyz_xxyz[k] = -g_x_0_xyyyz_xxyz[k] * cd_y[k] + g_x_0_xyyyz_xxyyz[k];

                g_x_0_xyyyyz_xxzz[k] = -g_x_0_xyyyz_xxzz[k] * cd_y[k] + g_x_0_xyyyz_xxyzz[k];

                g_x_0_xyyyyz_xyyy[k] = -g_x_0_xyyyz_xyyy[k] * cd_y[k] + g_x_0_xyyyz_xyyyy[k];

                g_x_0_xyyyyz_xyyz[k] = -g_x_0_xyyyz_xyyz[k] * cd_y[k] + g_x_0_xyyyz_xyyyz[k];

                g_x_0_xyyyyz_xyzz[k] = -g_x_0_xyyyz_xyzz[k] * cd_y[k] + g_x_0_xyyyz_xyyzz[k];

                g_x_0_xyyyyz_xzzz[k] = -g_x_0_xyyyz_xzzz[k] * cd_y[k] + g_x_0_xyyyz_xyzzz[k];

                g_x_0_xyyyyz_yyyy[k] = -g_x_0_xyyyz_yyyy[k] * cd_y[k] + g_x_0_xyyyz_yyyyy[k];

                g_x_0_xyyyyz_yyyz[k] = -g_x_0_xyyyz_yyyz[k] * cd_y[k] + g_x_0_xyyyz_yyyyz[k];

                g_x_0_xyyyyz_yyzz[k] = -g_x_0_xyyyz_yyzz[k] * cd_y[k] + g_x_0_xyyyz_yyyzz[k];

                g_x_0_xyyyyz_yzzz[k] = -g_x_0_xyyyz_yzzz[k] * cd_y[k] + g_x_0_xyyyz_yyzzz[k];

                g_x_0_xyyyyz_zzzz[k] = -g_x_0_xyyyz_zzzz[k] * cd_y[k] + g_x_0_xyyyz_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_xyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_xyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_xyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_xyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_xyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_xyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_xyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_xyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_xyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_xyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_xyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_xyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_xyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_xyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 269);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyzz_xxxx, g_x_0_xyyyzz_xxxy, g_x_0_xyyyzz_xxxz, g_x_0_xyyyzz_xxyy, g_x_0_xyyyzz_xxyz, g_x_0_xyyyzz_xxzz, g_x_0_xyyyzz_xyyy, g_x_0_xyyyzz_xyyz, g_x_0_xyyyzz_xyzz, g_x_0_xyyyzz_xzzz, g_x_0_xyyyzz_yyyy, g_x_0_xyyyzz_yyyz, g_x_0_xyyyzz_yyzz, g_x_0_xyyyzz_yzzz, g_x_0_xyyyzz_zzzz, g_x_0_xyyzz_xxxx, g_x_0_xyyzz_xxxxy, g_x_0_xyyzz_xxxy, g_x_0_xyyzz_xxxyy, g_x_0_xyyzz_xxxyz, g_x_0_xyyzz_xxxz, g_x_0_xyyzz_xxyy, g_x_0_xyyzz_xxyyy, g_x_0_xyyzz_xxyyz, g_x_0_xyyzz_xxyz, g_x_0_xyyzz_xxyzz, g_x_0_xyyzz_xxzz, g_x_0_xyyzz_xyyy, g_x_0_xyyzz_xyyyy, g_x_0_xyyzz_xyyyz, g_x_0_xyyzz_xyyz, g_x_0_xyyzz_xyyzz, g_x_0_xyyzz_xyzz, g_x_0_xyyzz_xyzzz, g_x_0_xyyzz_xzzz, g_x_0_xyyzz_yyyy, g_x_0_xyyzz_yyyyy, g_x_0_xyyzz_yyyyz, g_x_0_xyyzz_yyyz, g_x_0_xyyzz_yyyzz, g_x_0_xyyzz_yyzz, g_x_0_xyyzz_yyzzz, g_x_0_xyyzz_yzzz, g_x_0_xyyzz_yzzzz, g_x_0_xyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_xxxx[k] = -g_x_0_xyyzz_xxxx[k] * cd_y[k] + g_x_0_xyyzz_xxxxy[k];

                g_x_0_xyyyzz_xxxy[k] = -g_x_0_xyyzz_xxxy[k] * cd_y[k] + g_x_0_xyyzz_xxxyy[k];

                g_x_0_xyyyzz_xxxz[k] = -g_x_0_xyyzz_xxxz[k] * cd_y[k] + g_x_0_xyyzz_xxxyz[k];

                g_x_0_xyyyzz_xxyy[k] = -g_x_0_xyyzz_xxyy[k] * cd_y[k] + g_x_0_xyyzz_xxyyy[k];

                g_x_0_xyyyzz_xxyz[k] = -g_x_0_xyyzz_xxyz[k] * cd_y[k] + g_x_0_xyyzz_xxyyz[k];

                g_x_0_xyyyzz_xxzz[k] = -g_x_0_xyyzz_xxzz[k] * cd_y[k] + g_x_0_xyyzz_xxyzz[k];

                g_x_0_xyyyzz_xyyy[k] = -g_x_0_xyyzz_xyyy[k] * cd_y[k] + g_x_0_xyyzz_xyyyy[k];

                g_x_0_xyyyzz_xyyz[k] = -g_x_0_xyyzz_xyyz[k] * cd_y[k] + g_x_0_xyyzz_xyyyz[k];

                g_x_0_xyyyzz_xyzz[k] = -g_x_0_xyyzz_xyzz[k] * cd_y[k] + g_x_0_xyyzz_xyyzz[k];

                g_x_0_xyyyzz_xzzz[k] = -g_x_0_xyyzz_xzzz[k] * cd_y[k] + g_x_0_xyyzz_xyzzz[k];

                g_x_0_xyyyzz_yyyy[k] = -g_x_0_xyyzz_yyyy[k] * cd_y[k] + g_x_0_xyyzz_yyyyy[k];

                g_x_0_xyyyzz_yyyz[k] = -g_x_0_xyyzz_yyyz[k] * cd_y[k] + g_x_0_xyyzz_yyyyz[k];

                g_x_0_xyyyzz_yyzz[k] = -g_x_0_xyyzz_yyzz[k] * cd_y[k] + g_x_0_xyyzz_yyyzz[k];

                g_x_0_xyyyzz_yzzz[k] = -g_x_0_xyyzz_yzzz[k] * cd_y[k] + g_x_0_xyyzz_yyzzz[k];

                g_x_0_xyyyzz_zzzz[k] = -g_x_0_xyyzz_zzzz[k] * cd_y[k] + g_x_0_xyyzz_yzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_xyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_xyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_xyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_xyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_xyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_xyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_xyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_xyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_xyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_xyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_xyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_xyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_xyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_xyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 284);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzzz_xxxx, g_x_0_xyyzzz_xxxy, g_x_0_xyyzzz_xxxz, g_x_0_xyyzzz_xxyy, g_x_0_xyyzzz_xxyz, g_x_0_xyyzzz_xxzz, g_x_0_xyyzzz_xyyy, g_x_0_xyyzzz_xyyz, g_x_0_xyyzzz_xyzz, g_x_0_xyyzzz_xzzz, g_x_0_xyyzzz_yyyy, g_x_0_xyyzzz_yyyz, g_x_0_xyyzzz_yyzz, g_x_0_xyyzzz_yzzz, g_x_0_xyyzzz_zzzz, g_x_0_xyzzz_xxxx, g_x_0_xyzzz_xxxxy, g_x_0_xyzzz_xxxy, g_x_0_xyzzz_xxxyy, g_x_0_xyzzz_xxxyz, g_x_0_xyzzz_xxxz, g_x_0_xyzzz_xxyy, g_x_0_xyzzz_xxyyy, g_x_0_xyzzz_xxyyz, g_x_0_xyzzz_xxyz, g_x_0_xyzzz_xxyzz, g_x_0_xyzzz_xxzz, g_x_0_xyzzz_xyyy, g_x_0_xyzzz_xyyyy, g_x_0_xyzzz_xyyyz, g_x_0_xyzzz_xyyz, g_x_0_xyzzz_xyyzz, g_x_0_xyzzz_xyzz, g_x_0_xyzzz_xyzzz, g_x_0_xyzzz_xzzz, g_x_0_xyzzz_yyyy, g_x_0_xyzzz_yyyyy, g_x_0_xyzzz_yyyyz, g_x_0_xyzzz_yyyz, g_x_0_xyzzz_yyyzz, g_x_0_xyzzz_yyzz, g_x_0_xyzzz_yyzzz, g_x_0_xyzzz_yzzz, g_x_0_xyzzz_yzzzz, g_x_0_xyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_xxxx[k] = -g_x_0_xyzzz_xxxx[k] * cd_y[k] + g_x_0_xyzzz_xxxxy[k];

                g_x_0_xyyzzz_xxxy[k] = -g_x_0_xyzzz_xxxy[k] * cd_y[k] + g_x_0_xyzzz_xxxyy[k];

                g_x_0_xyyzzz_xxxz[k] = -g_x_0_xyzzz_xxxz[k] * cd_y[k] + g_x_0_xyzzz_xxxyz[k];

                g_x_0_xyyzzz_xxyy[k] = -g_x_0_xyzzz_xxyy[k] * cd_y[k] + g_x_0_xyzzz_xxyyy[k];

                g_x_0_xyyzzz_xxyz[k] = -g_x_0_xyzzz_xxyz[k] * cd_y[k] + g_x_0_xyzzz_xxyyz[k];

                g_x_0_xyyzzz_xxzz[k] = -g_x_0_xyzzz_xxzz[k] * cd_y[k] + g_x_0_xyzzz_xxyzz[k];

                g_x_0_xyyzzz_xyyy[k] = -g_x_0_xyzzz_xyyy[k] * cd_y[k] + g_x_0_xyzzz_xyyyy[k];

                g_x_0_xyyzzz_xyyz[k] = -g_x_0_xyzzz_xyyz[k] * cd_y[k] + g_x_0_xyzzz_xyyyz[k];

                g_x_0_xyyzzz_xyzz[k] = -g_x_0_xyzzz_xyzz[k] * cd_y[k] + g_x_0_xyzzz_xyyzz[k];

                g_x_0_xyyzzz_xzzz[k] = -g_x_0_xyzzz_xzzz[k] * cd_y[k] + g_x_0_xyzzz_xyzzz[k];

                g_x_0_xyyzzz_yyyy[k] = -g_x_0_xyzzz_yyyy[k] * cd_y[k] + g_x_0_xyzzz_yyyyy[k];

                g_x_0_xyyzzz_yyyz[k] = -g_x_0_xyzzz_yyyz[k] * cd_y[k] + g_x_0_xyzzz_yyyyz[k];

                g_x_0_xyyzzz_yyzz[k] = -g_x_0_xyzzz_yyzz[k] * cd_y[k] + g_x_0_xyzzz_yyyzz[k];

                g_x_0_xyyzzz_yzzz[k] = -g_x_0_xyzzz_yzzz[k] * cd_y[k] + g_x_0_xyzzz_yyzzz[k];

                g_x_0_xyyzzz_zzzz[k] = -g_x_0_xyzzz_zzzz[k] * cd_y[k] + g_x_0_xyzzz_yzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_xyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_xyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_xyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_xyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_xyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_xyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_xyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_xyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_xyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_xyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_xyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_xyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_xyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_xyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 299);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzzz_xxxx, g_x_0_xyzzzz_xxxy, g_x_0_xyzzzz_xxxz, g_x_0_xyzzzz_xxyy, g_x_0_xyzzzz_xxyz, g_x_0_xyzzzz_xxzz, g_x_0_xyzzzz_xyyy, g_x_0_xyzzzz_xyyz, g_x_0_xyzzzz_xyzz, g_x_0_xyzzzz_xzzz, g_x_0_xyzzzz_yyyy, g_x_0_xyzzzz_yyyz, g_x_0_xyzzzz_yyzz, g_x_0_xyzzzz_yzzz, g_x_0_xyzzzz_zzzz, g_x_0_xzzzz_xxxx, g_x_0_xzzzz_xxxxy, g_x_0_xzzzz_xxxy, g_x_0_xzzzz_xxxyy, g_x_0_xzzzz_xxxyz, g_x_0_xzzzz_xxxz, g_x_0_xzzzz_xxyy, g_x_0_xzzzz_xxyyy, g_x_0_xzzzz_xxyyz, g_x_0_xzzzz_xxyz, g_x_0_xzzzz_xxyzz, g_x_0_xzzzz_xxzz, g_x_0_xzzzz_xyyy, g_x_0_xzzzz_xyyyy, g_x_0_xzzzz_xyyyz, g_x_0_xzzzz_xyyz, g_x_0_xzzzz_xyyzz, g_x_0_xzzzz_xyzz, g_x_0_xzzzz_xyzzz, g_x_0_xzzzz_xzzz, g_x_0_xzzzz_yyyy, g_x_0_xzzzz_yyyyy, g_x_0_xzzzz_yyyyz, g_x_0_xzzzz_yyyz, g_x_0_xzzzz_yyyzz, g_x_0_xzzzz_yyzz, g_x_0_xzzzz_yyzzz, g_x_0_xzzzz_yzzz, g_x_0_xzzzz_yzzzz, g_x_0_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_xxxx[k] = -g_x_0_xzzzz_xxxx[k] * cd_y[k] + g_x_0_xzzzz_xxxxy[k];

                g_x_0_xyzzzz_xxxy[k] = -g_x_0_xzzzz_xxxy[k] * cd_y[k] + g_x_0_xzzzz_xxxyy[k];

                g_x_0_xyzzzz_xxxz[k] = -g_x_0_xzzzz_xxxz[k] * cd_y[k] + g_x_0_xzzzz_xxxyz[k];

                g_x_0_xyzzzz_xxyy[k] = -g_x_0_xzzzz_xxyy[k] * cd_y[k] + g_x_0_xzzzz_xxyyy[k];

                g_x_0_xyzzzz_xxyz[k] = -g_x_0_xzzzz_xxyz[k] * cd_y[k] + g_x_0_xzzzz_xxyyz[k];

                g_x_0_xyzzzz_xxzz[k] = -g_x_0_xzzzz_xxzz[k] * cd_y[k] + g_x_0_xzzzz_xxyzz[k];

                g_x_0_xyzzzz_xyyy[k] = -g_x_0_xzzzz_xyyy[k] * cd_y[k] + g_x_0_xzzzz_xyyyy[k];

                g_x_0_xyzzzz_xyyz[k] = -g_x_0_xzzzz_xyyz[k] * cd_y[k] + g_x_0_xzzzz_xyyyz[k];

                g_x_0_xyzzzz_xyzz[k] = -g_x_0_xzzzz_xyzz[k] * cd_y[k] + g_x_0_xzzzz_xyyzz[k];

                g_x_0_xyzzzz_xzzz[k] = -g_x_0_xzzzz_xzzz[k] * cd_y[k] + g_x_0_xzzzz_xyzzz[k];

                g_x_0_xyzzzz_yyyy[k] = -g_x_0_xzzzz_yyyy[k] * cd_y[k] + g_x_0_xzzzz_yyyyy[k];

                g_x_0_xyzzzz_yyyz[k] = -g_x_0_xzzzz_yyyz[k] * cd_y[k] + g_x_0_xzzzz_yyyyz[k];

                g_x_0_xyzzzz_yyzz[k] = -g_x_0_xzzzz_yyzz[k] * cd_y[k] + g_x_0_xzzzz_yyyzz[k];

                g_x_0_xyzzzz_yzzz[k] = -g_x_0_xzzzz_yzzz[k] * cd_y[k] + g_x_0_xzzzz_yyzzz[k];

                g_x_0_xyzzzz_zzzz[k] = -g_x_0_xzzzz_zzzz[k] * cd_y[k] + g_x_0_xzzzz_yzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_xzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_xzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_xzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_xzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_xzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_xzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_xzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_xzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_xzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_xzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_xzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_xzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_xzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_xzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_z, g_x_0_xzzzz_xxxx, g_x_0_xzzzz_xxxxz, g_x_0_xzzzz_xxxy, g_x_0_xzzzz_xxxyz, g_x_0_xzzzz_xxxz, g_x_0_xzzzz_xxxzz, g_x_0_xzzzz_xxyy, g_x_0_xzzzz_xxyyz, g_x_0_xzzzz_xxyz, g_x_0_xzzzz_xxyzz, g_x_0_xzzzz_xxzz, g_x_0_xzzzz_xxzzz, g_x_0_xzzzz_xyyy, g_x_0_xzzzz_xyyyz, g_x_0_xzzzz_xyyz, g_x_0_xzzzz_xyyzz, g_x_0_xzzzz_xyzz, g_x_0_xzzzz_xyzzz, g_x_0_xzzzz_xzzz, g_x_0_xzzzz_xzzzz, g_x_0_xzzzz_yyyy, g_x_0_xzzzz_yyyyz, g_x_0_xzzzz_yyyz, g_x_0_xzzzz_yyyzz, g_x_0_xzzzz_yyzz, g_x_0_xzzzz_yyzzz, g_x_0_xzzzz_yzzz, g_x_0_xzzzz_yzzzz, g_x_0_xzzzz_zzzz, g_x_0_xzzzz_zzzzz, g_x_0_xzzzzz_xxxx, g_x_0_xzzzzz_xxxy, g_x_0_xzzzzz_xxxz, g_x_0_xzzzzz_xxyy, g_x_0_xzzzzz_xxyz, g_x_0_xzzzzz_xxzz, g_x_0_xzzzzz_xyyy, g_x_0_xzzzzz_xyyz, g_x_0_xzzzzz_xyzz, g_x_0_xzzzzz_xzzz, g_x_0_xzzzzz_yyyy, g_x_0_xzzzzz_yyyz, g_x_0_xzzzzz_yyzz, g_x_0_xzzzzz_yzzz, g_x_0_xzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_xxxx[k] = -g_x_0_xzzzz_xxxx[k] * cd_z[k] + g_x_0_xzzzz_xxxxz[k];

                g_x_0_xzzzzz_xxxy[k] = -g_x_0_xzzzz_xxxy[k] * cd_z[k] + g_x_0_xzzzz_xxxyz[k];

                g_x_0_xzzzzz_xxxz[k] = -g_x_0_xzzzz_xxxz[k] * cd_z[k] + g_x_0_xzzzz_xxxzz[k];

                g_x_0_xzzzzz_xxyy[k] = -g_x_0_xzzzz_xxyy[k] * cd_z[k] + g_x_0_xzzzz_xxyyz[k];

                g_x_0_xzzzzz_xxyz[k] = -g_x_0_xzzzz_xxyz[k] * cd_z[k] + g_x_0_xzzzz_xxyzz[k];

                g_x_0_xzzzzz_xxzz[k] = -g_x_0_xzzzz_xxzz[k] * cd_z[k] + g_x_0_xzzzz_xxzzz[k];

                g_x_0_xzzzzz_xyyy[k] = -g_x_0_xzzzz_xyyy[k] * cd_z[k] + g_x_0_xzzzz_xyyyz[k];

                g_x_0_xzzzzz_xyyz[k] = -g_x_0_xzzzz_xyyz[k] * cd_z[k] + g_x_0_xzzzz_xyyzz[k];

                g_x_0_xzzzzz_xyzz[k] = -g_x_0_xzzzz_xyzz[k] * cd_z[k] + g_x_0_xzzzz_xyzzz[k];

                g_x_0_xzzzzz_xzzz[k] = -g_x_0_xzzzz_xzzz[k] * cd_z[k] + g_x_0_xzzzz_xzzzz[k];

                g_x_0_xzzzzz_yyyy[k] = -g_x_0_xzzzz_yyyy[k] * cd_z[k] + g_x_0_xzzzz_yyyyz[k];

                g_x_0_xzzzzz_yyyz[k] = -g_x_0_xzzzz_yyyz[k] * cd_z[k] + g_x_0_xzzzz_yyyzz[k];

                g_x_0_xzzzzz_yyzz[k] = -g_x_0_xzzzz_yyzz[k] * cd_z[k] + g_x_0_xzzzz_yyzzz[k];

                g_x_0_xzzzzz_yzzz[k] = -g_x_0_xzzzz_yzzz[k] * cd_z[k] + g_x_0_xzzzz_yzzzz[k];

                g_x_0_xzzzzz_zzzz[k] = -g_x_0_xzzzz_zzzz[k] * cd_z[k] + g_x_0_xzzzz_zzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_yyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_yyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_yyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_yyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_yyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_yyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_yyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_yyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_yyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_yyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_yyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_yyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_yyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_yyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 329);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyy_xxxx, g_x_0_yyyyy_xxxxy, g_x_0_yyyyy_xxxy, g_x_0_yyyyy_xxxyy, g_x_0_yyyyy_xxxyz, g_x_0_yyyyy_xxxz, g_x_0_yyyyy_xxyy, g_x_0_yyyyy_xxyyy, g_x_0_yyyyy_xxyyz, g_x_0_yyyyy_xxyz, g_x_0_yyyyy_xxyzz, g_x_0_yyyyy_xxzz, g_x_0_yyyyy_xyyy, g_x_0_yyyyy_xyyyy, g_x_0_yyyyy_xyyyz, g_x_0_yyyyy_xyyz, g_x_0_yyyyy_xyyzz, g_x_0_yyyyy_xyzz, g_x_0_yyyyy_xyzzz, g_x_0_yyyyy_xzzz, g_x_0_yyyyy_yyyy, g_x_0_yyyyy_yyyyy, g_x_0_yyyyy_yyyyz, g_x_0_yyyyy_yyyz, g_x_0_yyyyy_yyyzz, g_x_0_yyyyy_yyzz, g_x_0_yyyyy_yyzzz, g_x_0_yyyyy_yzzz, g_x_0_yyyyy_yzzzz, g_x_0_yyyyy_zzzz, g_x_0_yyyyyy_xxxx, g_x_0_yyyyyy_xxxy, g_x_0_yyyyyy_xxxz, g_x_0_yyyyyy_xxyy, g_x_0_yyyyyy_xxyz, g_x_0_yyyyyy_xxzz, g_x_0_yyyyyy_xyyy, g_x_0_yyyyyy_xyyz, g_x_0_yyyyyy_xyzz, g_x_0_yyyyyy_xzzz, g_x_0_yyyyyy_yyyy, g_x_0_yyyyyy_yyyz, g_x_0_yyyyyy_yyzz, g_x_0_yyyyyy_yzzz, g_x_0_yyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_xxxx[k] = -g_x_0_yyyyy_xxxx[k] * cd_y[k] + g_x_0_yyyyy_xxxxy[k];

                g_x_0_yyyyyy_xxxy[k] = -g_x_0_yyyyy_xxxy[k] * cd_y[k] + g_x_0_yyyyy_xxxyy[k];

                g_x_0_yyyyyy_xxxz[k] = -g_x_0_yyyyy_xxxz[k] * cd_y[k] + g_x_0_yyyyy_xxxyz[k];

                g_x_0_yyyyyy_xxyy[k] = -g_x_0_yyyyy_xxyy[k] * cd_y[k] + g_x_0_yyyyy_xxyyy[k];

                g_x_0_yyyyyy_xxyz[k] = -g_x_0_yyyyy_xxyz[k] * cd_y[k] + g_x_0_yyyyy_xxyyz[k];

                g_x_0_yyyyyy_xxzz[k] = -g_x_0_yyyyy_xxzz[k] * cd_y[k] + g_x_0_yyyyy_xxyzz[k];

                g_x_0_yyyyyy_xyyy[k] = -g_x_0_yyyyy_xyyy[k] * cd_y[k] + g_x_0_yyyyy_xyyyy[k];

                g_x_0_yyyyyy_xyyz[k] = -g_x_0_yyyyy_xyyz[k] * cd_y[k] + g_x_0_yyyyy_xyyyz[k];

                g_x_0_yyyyyy_xyzz[k] = -g_x_0_yyyyy_xyzz[k] * cd_y[k] + g_x_0_yyyyy_xyyzz[k];

                g_x_0_yyyyyy_xzzz[k] = -g_x_0_yyyyy_xzzz[k] * cd_y[k] + g_x_0_yyyyy_xyzzz[k];

                g_x_0_yyyyyy_yyyy[k] = -g_x_0_yyyyy_yyyy[k] * cd_y[k] + g_x_0_yyyyy_yyyyy[k];

                g_x_0_yyyyyy_yyyz[k] = -g_x_0_yyyyy_yyyz[k] * cd_y[k] + g_x_0_yyyyy_yyyyz[k];

                g_x_0_yyyyyy_yyzz[k] = -g_x_0_yyyyy_yyzz[k] * cd_y[k] + g_x_0_yyyyy_yyyzz[k];

                g_x_0_yyyyyy_yzzz[k] = -g_x_0_yyyyy_yzzz[k] * cd_y[k] + g_x_0_yyyyy_yyzzz[k];

                g_x_0_yyyyyy_zzzz[k] = -g_x_0_yyyyy_zzzz[k] * cd_y[k] + g_x_0_yyyyy_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_yyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_yyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_yyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_yyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_yyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 335);

            auto g_x_0_yyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_yyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_yyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_yyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_yyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_yyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_yyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_yyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_yyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 344);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyyz_xxxx, g_x_0_yyyyyz_xxxy, g_x_0_yyyyyz_xxxz, g_x_0_yyyyyz_xxyy, g_x_0_yyyyyz_xxyz, g_x_0_yyyyyz_xxzz, g_x_0_yyyyyz_xyyy, g_x_0_yyyyyz_xyyz, g_x_0_yyyyyz_xyzz, g_x_0_yyyyyz_xzzz, g_x_0_yyyyyz_yyyy, g_x_0_yyyyyz_yyyz, g_x_0_yyyyyz_yyzz, g_x_0_yyyyyz_yzzz, g_x_0_yyyyyz_zzzz, g_x_0_yyyyz_xxxx, g_x_0_yyyyz_xxxxy, g_x_0_yyyyz_xxxy, g_x_0_yyyyz_xxxyy, g_x_0_yyyyz_xxxyz, g_x_0_yyyyz_xxxz, g_x_0_yyyyz_xxyy, g_x_0_yyyyz_xxyyy, g_x_0_yyyyz_xxyyz, g_x_0_yyyyz_xxyz, g_x_0_yyyyz_xxyzz, g_x_0_yyyyz_xxzz, g_x_0_yyyyz_xyyy, g_x_0_yyyyz_xyyyy, g_x_0_yyyyz_xyyyz, g_x_0_yyyyz_xyyz, g_x_0_yyyyz_xyyzz, g_x_0_yyyyz_xyzz, g_x_0_yyyyz_xyzzz, g_x_0_yyyyz_xzzz, g_x_0_yyyyz_yyyy, g_x_0_yyyyz_yyyyy, g_x_0_yyyyz_yyyyz, g_x_0_yyyyz_yyyz, g_x_0_yyyyz_yyyzz, g_x_0_yyyyz_yyzz, g_x_0_yyyyz_yyzzz, g_x_0_yyyyz_yzzz, g_x_0_yyyyz_yzzzz, g_x_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_xxxx[k] = -g_x_0_yyyyz_xxxx[k] * cd_y[k] + g_x_0_yyyyz_xxxxy[k];

                g_x_0_yyyyyz_xxxy[k] = -g_x_0_yyyyz_xxxy[k] * cd_y[k] + g_x_0_yyyyz_xxxyy[k];

                g_x_0_yyyyyz_xxxz[k] = -g_x_0_yyyyz_xxxz[k] * cd_y[k] + g_x_0_yyyyz_xxxyz[k];

                g_x_0_yyyyyz_xxyy[k] = -g_x_0_yyyyz_xxyy[k] * cd_y[k] + g_x_0_yyyyz_xxyyy[k];

                g_x_0_yyyyyz_xxyz[k] = -g_x_0_yyyyz_xxyz[k] * cd_y[k] + g_x_0_yyyyz_xxyyz[k];

                g_x_0_yyyyyz_xxzz[k] = -g_x_0_yyyyz_xxzz[k] * cd_y[k] + g_x_0_yyyyz_xxyzz[k];

                g_x_0_yyyyyz_xyyy[k] = -g_x_0_yyyyz_xyyy[k] * cd_y[k] + g_x_0_yyyyz_xyyyy[k];

                g_x_0_yyyyyz_xyyz[k] = -g_x_0_yyyyz_xyyz[k] * cd_y[k] + g_x_0_yyyyz_xyyyz[k];

                g_x_0_yyyyyz_xyzz[k] = -g_x_0_yyyyz_xyzz[k] * cd_y[k] + g_x_0_yyyyz_xyyzz[k];

                g_x_0_yyyyyz_xzzz[k] = -g_x_0_yyyyz_xzzz[k] * cd_y[k] + g_x_0_yyyyz_xyzzz[k];

                g_x_0_yyyyyz_yyyy[k] = -g_x_0_yyyyz_yyyy[k] * cd_y[k] + g_x_0_yyyyz_yyyyy[k];

                g_x_0_yyyyyz_yyyz[k] = -g_x_0_yyyyz_yyyz[k] * cd_y[k] + g_x_0_yyyyz_yyyyz[k];

                g_x_0_yyyyyz_yyzz[k] = -g_x_0_yyyyz_yyzz[k] * cd_y[k] + g_x_0_yyyyz_yyyzz[k];

                g_x_0_yyyyyz_yzzz[k] = -g_x_0_yyyyz_yzzz[k] * cd_y[k] + g_x_0_yyyyz_yyzzz[k];

                g_x_0_yyyyyz_zzzz[k] = -g_x_0_yyyyz_zzzz[k] * cd_y[k] + g_x_0_yyyyz_yzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_yyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_yyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_yyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_yyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_yyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_yyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_yyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_yyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_yyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_yyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_yyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 356);

            auto g_x_0_yyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_yyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_yyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 359);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyzz_xxxx, g_x_0_yyyyzz_xxxy, g_x_0_yyyyzz_xxxz, g_x_0_yyyyzz_xxyy, g_x_0_yyyyzz_xxyz, g_x_0_yyyyzz_xxzz, g_x_0_yyyyzz_xyyy, g_x_0_yyyyzz_xyyz, g_x_0_yyyyzz_xyzz, g_x_0_yyyyzz_xzzz, g_x_0_yyyyzz_yyyy, g_x_0_yyyyzz_yyyz, g_x_0_yyyyzz_yyzz, g_x_0_yyyyzz_yzzz, g_x_0_yyyyzz_zzzz, g_x_0_yyyzz_xxxx, g_x_0_yyyzz_xxxxy, g_x_0_yyyzz_xxxy, g_x_0_yyyzz_xxxyy, g_x_0_yyyzz_xxxyz, g_x_0_yyyzz_xxxz, g_x_0_yyyzz_xxyy, g_x_0_yyyzz_xxyyy, g_x_0_yyyzz_xxyyz, g_x_0_yyyzz_xxyz, g_x_0_yyyzz_xxyzz, g_x_0_yyyzz_xxzz, g_x_0_yyyzz_xyyy, g_x_0_yyyzz_xyyyy, g_x_0_yyyzz_xyyyz, g_x_0_yyyzz_xyyz, g_x_0_yyyzz_xyyzz, g_x_0_yyyzz_xyzz, g_x_0_yyyzz_xyzzz, g_x_0_yyyzz_xzzz, g_x_0_yyyzz_yyyy, g_x_0_yyyzz_yyyyy, g_x_0_yyyzz_yyyyz, g_x_0_yyyzz_yyyz, g_x_0_yyyzz_yyyzz, g_x_0_yyyzz_yyzz, g_x_0_yyyzz_yyzzz, g_x_0_yyyzz_yzzz, g_x_0_yyyzz_yzzzz, g_x_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_xxxx[k] = -g_x_0_yyyzz_xxxx[k] * cd_y[k] + g_x_0_yyyzz_xxxxy[k];

                g_x_0_yyyyzz_xxxy[k] = -g_x_0_yyyzz_xxxy[k] * cd_y[k] + g_x_0_yyyzz_xxxyy[k];

                g_x_0_yyyyzz_xxxz[k] = -g_x_0_yyyzz_xxxz[k] * cd_y[k] + g_x_0_yyyzz_xxxyz[k];

                g_x_0_yyyyzz_xxyy[k] = -g_x_0_yyyzz_xxyy[k] * cd_y[k] + g_x_0_yyyzz_xxyyy[k];

                g_x_0_yyyyzz_xxyz[k] = -g_x_0_yyyzz_xxyz[k] * cd_y[k] + g_x_0_yyyzz_xxyyz[k];

                g_x_0_yyyyzz_xxzz[k] = -g_x_0_yyyzz_xxzz[k] * cd_y[k] + g_x_0_yyyzz_xxyzz[k];

                g_x_0_yyyyzz_xyyy[k] = -g_x_0_yyyzz_xyyy[k] * cd_y[k] + g_x_0_yyyzz_xyyyy[k];

                g_x_0_yyyyzz_xyyz[k] = -g_x_0_yyyzz_xyyz[k] * cd_y[k] + g_x_0_yyyzz_xyyyz[k];

                g_x_0_yyyyzz_xyzz[k] = -g_x_0_yyyzz_xyzz[k] * cd_y[k] + g_x_0_yyyzz_xyyzz[k];

                g_x_0_yyyyzz_xzzz[k] = -g_x_0_yyyzz_xzzz[k] * cd_y[k] + g_x_0_yyyzz_xyzzz[k];

                g_x_0_yyyyzz_yyyy[k] = -g_x_0_yyyzz_yyyy[k] * cd_y[k] + g_x_0_yyyzz_yyyyy[k];

                g_x_0_yyyyzz_yyyz[k] = -g_x_0_yyyzz_yyyz[k] * cd_y[k] + g_x_0_yyyzz_yyyyz[k];

                g_x_0_yyyyzz_yyzz[k] = -g_x_0_yyyzz_yyzz[k] * cd_y[k] + g_x_0_yyyzz_yyyzz[k];

                g_x_0_yyyyzz_yzzz[k] = -g_x_0_yyyzz_yzzz[k] * cd_y[k] + g_x_0_yyyzz_yyzzz[k];

                g_x_0_yyyyzz_zzzz[k] = -g_x_0_yyyzz_zzzz[k] * cd_y[k] + g_x_0_yyyzz_yzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 360);

            auto g_x_0_yyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 361);

            auto g_x_0_yyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 362);

            auto g_x_0_yyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 363);

            auto g_x_0_yyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 364);

            auto g_x_0_yyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 365);

            auto g_x_0_yyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 366);

            auto g_x_0_yyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 367);

            auto g_x_0_yyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 368);

            auto g_x_0_yyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 369);

            auto g_x_0_yyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 370);

            auto g_x_0_yyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 371);

            auto g_x_0_yyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 372);

            auto g_x_0_yyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 373);

            auto g_x_0_yyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 374);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzzz_xxxx, g_x_0_yyyzzz_xxxy, g_x_0_yyyzzz_xxxz, g_x_0_yyyzzz_xxyy, g_x_0_yyyzzz_xxyz, g_x_0_yyyzzz_xxzz, g_x_0_yyyzzz_xyyy, g_x_0_yyyzzz_xyyz, g_x_0_yyyzzz_xyzz, g_x_0_yyyzzz_xzzz, g_x_0_yyyzzz_yyyy, g_x_0_yyyzzz_yyyz, g_x_0_yyyzzz_yyzz, g_x_0_yyyzzz_yzzz, g_x_0_yyyzzz_zzzz, g_x_0_yyzzz_xxxx, g_x_0_yyzzz_xxxxy, g_x_0_yyzzz_xxxy, g_x_0_yyzzz_xxxyy, g_x_0_yyzzz_xxxyz, g_x_0_yyzzz_xxxz, g_x_0_yyzzz_xxyy, g_x_0_yyzzz_xxyyy, g_x_0_yyzzz_xxyyz, g_x_0_yyzzz_xxyz, g_x_0_yyzzz_xxyzz, g_x_0_yyzzz_xxzz, g_x_0_yyzzz_xyyy, g_x_0_yyzzz_xyyyy, g_x_0_yyzzz_xyyyz, g_x_0_yyzzz_xyyz, g_x_0_yyzzz_xyyzz, g_x_0_yyzzz_xyzz, g_x_0_yyzzz_xyzzz, g_x_0_yyzzz_xzzz, g_x_0_yyzzz_yyyy, g_x_0_yyzzz_yyyyy, g_x_0_yyzzz_yyyyz, g_x_0_yyzzz_yyyz, g_x_0_yyzzz_yyyzz, g_x_0_yyzzz_yyzz, g_x_0_yyzzz_yyzzz, g_x_0_yyzzz_yzzz, g_x_0_yyzzz_yzzzz, g_x_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_xxxx[k] = -g_x_0_yyzzz_xxxx[k] * cd_y[k] + g_x_0_yyzzz_xxxxy[k];

                g_x_0_yyyzzz_xxxy[k] = -g_x_0_yyzzz_xxxy[k] * cd_y[k] + g_x_0_yyzzz_xxxyy[k];

                g_x_0_yyyzzz_xxxz[k] = -g_x_0_yyzzz_xxxz[k] * cd_y[k] + g_x_0_yyzzz_xxxyz[k];

                g_x_0_yyyzzz_xxyy[k] = -g_x_0_yyzzz_xxyy[k] * cd_y[k] + g_x_0_yyzzz_xxyyy[k];

                g_x_0_yyyzzz_xxyz[k] = -g_x_0_yyzzz_xxyz[k] * cd_y[k] + g_x_0_yyzzz_xxyyz[k];

                g_x_0_yyyzzz_xxzz[k] = -g_x_0_yyzzz_xxzz[k] * cd_y[k] + g_x_0_yyzzz_xxyzz[k];

                g_x_0_yyyzzz_xyyy[k] = -g_x_0_yyzzz_xyyy[k] * cd_y[k] + g_x_0_yyzzz_xyyyy[k];

                g_x_0_yyyzzz_xyyz[k] = -g_x_0_yyzzz_xyyz[k] * cd_y[k] + g_x_0_yyzzz_xyyyz[k];

                g_x_0_yyyzzz_xyzz[k] = -g_x_0_yyzzz_xyzz[k] * cd_y[k] + g_x_0_yyzzz_xyyzz[k];

                g_x_0_yyyzzz_xzzz[k] = -g_x_0_yyzzz_xzzz[k] * cd_y[k] + g_x_0_yyzzz_xyzzz[k];

                g_x_0_yyyzzz_yyyy[k] = -g_x_0_yyzzz_yyyy[k] * cd_y[k] + g_x_0_yyzzz_yyyyy[k];

                g_x_0_yyyzzz_yyyz[k] = -g_x_0_yyzzz_yyyz[k] * cd_y[k] + g_x_0_yyzzz_yyyyz[k];

                g_x_0_yyyzzz_yyzz[k] = -g_x_0_yyzzz_yyzz[k] * cd_y[k] + g_x_0_yyzzz_yyyzz[k];

                g_x_0_yyyzzz_yzzz[k] = -g_x_0_yyzzz_yzzz[k] * cd_y[k] + g_x_0_yyzzz_yyzzz[k];

                g_x_0_yyyzzz_zzzz[k] = -g_x_0_yyzzz_zzzz[k] * cd_y[k] + g_x_0_yyzzz_yzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 375);

            auto g_x_0_yyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 376);

            auto g_x_0_yyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 377);

            auto g_x_0_yyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 378);

            auto g_x_0_yyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 379);

            auto g_x_0_yyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 380);

            auto g_x_0_yyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 381);

            auto g_x_0_yyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 382);

            auto g_x_0_yyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 383);

            auto g_x_0_yyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 384);

            auto g_x_0_yyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 385);

            auto g_x_0_yyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 386);

            auto g_x_0_yyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 387);

            auto g_x_0_yyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 388);

            auto g_x_0_yyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 389);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzzz_xxxx, g_x_0_yyzzzz_xxxy, g_x_0_yyzzzz_xxxz, g_x_0_yyzzzz_xxyy, g_x_0_yyzzzz_xxyz, g_x_0_yyzzzz_xxzz, g_x_0_yyzzzz_xyyy, g_x_0_yyzzzz_xyyz, g_x_0_yyzzzz_xyzz, g_x_0_yyzzzz_xzzz, g_x_0_yyzzzz_yyyy, g_x_0_yyzzzz_yyyz, g_x_0_yyzzzz_yyzz, g_x_0_yyzzzz_yzzz, g_x_0_yyzzzz_zzzz, g_x_0_yzzzz_xxxx, g_x_0_yzzzz_xxxxy, g_x_0_yzzzz_xxxy, g_x_0_yzzzz_xxxyy, g_x_0_yzzzz_xxxyz, g_x_0_yzzzz_xxxz, g_x_0_yzzzz_xxyy, g_x_0_yzzzz_xxyyy, g_x_0_yzzzz_xxyyz, g_x_0_yzzzz_xxyz, g_x_0_yzzzz_xxyzz, g_x_0_yzzzz_xxzz, g_x_0_yzzzz_xyyy, g_x_0_yzzzz_xyyyy, g_x_0_yzzzz_xyyyz, g_x_0_yzzzz_xyyz, g_x_0_yzzzz_xyyzz, g_x_0_yzzzz_xyzz, g_x_0_yzzzz_xyzzz, g_x_0_yzzzz_xzzz, g_x_0_yzzzz_yyyy, g_x_0_yzzzz_yyyyy, g_x_0_yzzzz_yyyyz, g_x_0_yzzzz_yyyz, g_x_0_yzzzz_yyyzz, g_x_0_yzzzz_yyzz, g_x_0_yzzzz_yyzzz, g_x_0_yzzzz_yzzz, g_x_0_yzzzz_yzzzz, g_x_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_xxxx[k] = -g_x_0_yzzzz_xxxx[k] * cd_y[k] + g_x_0_yzzzz_xxxxy[k];

                g_x_0_yyzzzz_xxxy[k] = -g_x_0_yzzzz_xxxy[k] * cd_y[k] + g_x_0_yzzzz_xxxyy[k];

                g_x_0_yyzzzz_xxxz[k] = -g_x_0_yzzzz_xxxz[k] * cd_y[k] + g_x_0_yzzzz_xxxyz[k];

                g_x_0_yyzzzz_xxyy[k] = -g_x_0_yzzzz_xxyy[k] * cd_y[k] + g_x_0_yzzzz_xxyyy[k];

                g_x_0_yyzzzz_xxyz[k] = -g_x_0_yzzzz_xxyz[k] * cd_y[k] + g_x_0_yzzzz_xxyyz[k];

                g_x_0_yyzzzz_xxzz[k] = -g_x_0_yzzzz_xxzz[k] * cd_y[k] + g_x_0_yzzzz_xxyzz[k];

                g_x_0_yyzzzz_xyyy[k] = -g_x_0_yzzzz_xyyy[k] * cd_y[k] + g_x_0_yzzzz_xyyyy[k];

                g_x_0_yyzzzz_xyyz[k] = -g_x_0_yzzzz_xyyz[k] * cd_y[k] + g_x_0_yzzzz_xyyyz[k];

                g_x_0_yyzzzz_xyzz[k] = -g_x_0_yzzzz_xyzz[k] * cd_y[k] + g_x_0_yzzzz_xyyzz[k];

                g_x_0_yyzzzz_xzzz[k] = -g_x_0_yzzzz_xzzz[k] * cd_y[k] + g_x_0_yzzzz_xyzzz[k];

                g_x_0_yyzzzz_yyyy[k] = -g_x_0_yzzzz_yyyy[k] * cd_y[k] + g_x_0_yzzzz_yyyyy[k];

                g_x_0_yyzzzz_yyyz[k] = -g_x_0_yzzzz_yyyz[k] * cd_y[k] + g_x_0_yzzzz_yyyyz[k];

                g_x_0_yyzzzz_yyzz[k] = -g_x_0_yzzzz_yyzz[k] * cd_y[k] + g_x_0_yzzzz_yyyzz[k];

                g_x_0_yyzzzz_yzzz[k] = -g_x_0_yzzzz_yzzz[k] * cd_y[k] + g_x_0_yzzzz_yyzzz[k];

                g_x_0_yyzzzz_zzzz[k] = -g_x_0_yzzzz_zzzz[k] * cd_y[k] + g_x_0_yzzzz_yzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 390);

            auto g_x_0_yzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 391);

            auto g_x_0_yzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 392);

            auto g_x_0_yzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 393);

            auto g_x_0_yzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 394);

            auto g_x_0_yzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 395);

            auto g_x_0_yzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 396);

            auto g_x_0_yzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 397);

            auto g_x_0_yzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 398);

            auto g_x_0_yzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 399);

            auto g_x_0_yzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 400);

            auto g_x_0_yzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 401);

            auto g_x_0_yzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 402);

            auto g_x_0_yzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 403);

            auto g_x_0_yzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 404);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzzz_xxxx, g_x_0_yzzzzz_xxxy, g_x_0_yzzzzz_xxxz, g_x_0_yzzzzz_xxyy, g_x_0_yzzzzz_xxyz, g_x_0_yzzzzz_xxzz, g_x_0_yzzzzz_xyyy, g_x_0_yzzzzz_xyyz, g_x_0_yzzzzz_xyzz, g_x_0_yzzzzz_xzzz, g_x_0_yzzzzz_yyyy, g_x_0_yzzzzz_yyyz, g_x_0_yzzzzz_yyzz, g_x_0_yzzzzz_yzzz, g_x_0_yzzzzz_zzzz, g_x_0_zzzzz_xxxx, g_x_0_zzzzz_xxxxy, g_x_0_zzzzz_xxxy, g_x_0_zzzzz_xxxyy, g_x_0_zzzzz_xxxyz, g_x_0_zzzzz_xxxz, g_x_0_zzzzz_xxyy, g_x_0_zzzzz_xxyyy, g_x_0_zzzzz_xxyyz, g_x_0_zzzzz_xxyz, g_x_0_zzzzz_xxyzz, g_x_0_zzzzz_xxzz, g_x_0_zzzzz_xyyy, g_x_0_zzzzz_xyyyy, g_x_0_zzzzz_xyyyz, g_x_0_zzzzz_xyyz, g_x_0_zzzzz_xyyzz, g_x_0_zzzzz_xyzz, g_x_0_zzzzz_xyzzz, g_x_0_zzzzz_xzzz, g_x_0_zzzzz_yyyy, g_x_0_zzzzz_yyyyy, g_x_0_zzzzz_yyyyz, g_x_0_zzzzz_yyyz, g_x_0_zzzzz_yyyzz, g_x_0_zzzzz_yyzz, g_x_0_zzzzz_yyzzz, g_x_0_zzzzz_yzzz, g_x_0_zzzzz_yzzzz, g_x_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_xxxx[k] = -g_x_0_zzzzz_xxxx[k] * cd_y[k] + g_x_0_zzzzz_xxxxy[k];

                g_x_0_yzzzzz_xxxy[k] = -g_x_0_zzzzz_xxxy[k] * cd_y[k] + g_x_0_zzzzz_xxxyy[k];

                g_x_0_yzzzzz_xxxz[k] = -g_x_0_zzzzz_xxxz[k] * cd_y[k] + g_x_0_zzzzz_xxxyz[k];

                g_x_0_yzzzzz_xxyy[k] = -g_x_0_zzzzz_xxyy[k] * cd_y[k] + g_x_0_zzzzz_xxyyy[k];

                g_x_0_yzzzzz_xxyz[k] = -g_x_0_zzzzz_xxyz[k] * cd_y[k] + g_x_0_zzzzz_xxyyz[k];

                g_x_0_yzzzzz_xxzz[k] = -g_x_0_zzzzz_xxzz[k] * cd_y[k] + g_x_0_zzzzz_xxyzz[k];

                g_x_0_yzzzzz_xyyy[k] = -g_x_0_zzzzz_xyyy[k] * cd_y[k] + g_x_0_zzzzz_xyyyy[k];

                g_x_0_yzzzzz_xyyz[k] = -g_x_0_zzzzz_xyyz[k] * cd_y[k] + g_x_0_zzzzz_xyyyz[k];

                g_x_0_yzzzzz_xyzz[k] = -g_x_0_zzzzz_xyzz[k] * cd_y[k] + g_x_0_zzzzz_xyyzz[k];

                g_x_0_yzzzzz_xzzz[k] = -g_x_0_zzzzz_xzzz[k] * cd_y[k] + g_x_0_zzzzz_xyzzz[k];

                g_x_0_yzzzzz_yyyy[k] = -g_x_0_zzzzz_yyyy[k] * cd_y[k] + g_x_0_zzzzz_yyyyy[k];

                g_x_0_yzzzzz_yyyz[k] = -g_x_0_zzzzz_yyyz[k] * cd_y[k] + g_x_0_zzzzz_yyyyz[k];

                g_x_0_yzzzzz_yyzz[k] = -g_x_0_zzzzz_yyzz[k] * cd_y[k] + g_x_0_zzzzz_yyyzz[k];

                g_x_0_yzzzzz_yzzz[k] = -g_x_0_zzzzz_yzzz[k] * cd_y[k] + g_x_0_zzzzz_yyzzz[k];

                g_x_0_yzzzzz_zzzz[k] = -g_x_0_zzzzz_zzzz[k] * cd_y[k] + g_x_0_zzzzz_yzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 405);

            auto g_x_0_zzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 406);

            auto g_x_0_zzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 407);

            auto g_x_0_zzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 408);

            auto g_x_0_zzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 409);

            auto g_x_0_zzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 410);

            auto g_x_0_zzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 411);

            auto g_x_0_zzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 412);

            auto g_x_0_zzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 413);

            auto g_x_0_zzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 414);

            auto g_x_0_zzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 415);

            auto g_x_0_zzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 416);

            auto g_x_0_zzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 417);

            auto g_x_0_zzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 418);

            auto g_x_0_zzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 0 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_z, g_x_0_zzzzz_xxxx, g_x_0_zzzzz_xxxxz, g_x_0_zzzzz_xxxy, g_x_0_zzzzz_xxxyz, g_x_0_zzzzz_xxxz, g_x_0_zzzzz_xxxzz, g_x_0_zzzzz_xxyy, g_x_0_zzzzz_xxyyz, g_x_0_zzzzz_xxyz, g_x_0_zzzzz_xxyzz, g_x_0_zzzzz_xxzz, g_x_0_zzzzz_xxzzz, g_x_0_zzzzz_xyyy, g_x_0_zzzzz_xyyyz, g_x_0_zzzzz_xyyz, g_x_0_zzzzz_xyyzz, g_x_0_zzzzz_xyzz, g_x_0_zzzzz_xyzzz, g_x_0_zzzzz_xzzz, g_x_0_zzzzz_xzzzz, g_x_0_zzzzz_yyyy, g_x_0_zzzzz_yyyyz, g_x_0_zzzzz_yyyz, g_x_0_zzzzz_yyyzz, g_x_0_zzzzz_yyzz, g_x_0_zzzzz_yyzzz, g_x_0_zzzzz_yzzz, g_x_0_zzzzz_yzzzz, g_x_0_zzzzz_zzzz, g_x_0_zzzzz_zzzzz, g_x_0_zzzzzz_xxxx, g_x_0_zzzzzz_xxxy, g_x_0_zzzzzz_xxxz, g_x_0_zzzzzz_xxyy, g_x_0_zzzzzz_xxyz, g_x_0_zzzzzz_xxzz, g_x_0_zzzzzz_xyyy, g_x_0_zzzzzz_xyyz, g_x_0_zzzzzz_xyzz, g_x_0_zzzzzz_xzzz, g_x_0_zzzzzz_yyyy, g_x_0_zzzzzz_yyyz, g_x_0_zzzzzz_yyzz, g_x_0_zzzzzz_yzzz, g_x_0_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_xxxx[k] = -g_x_0_zzzzz_xxxx[k] * cd_z[k] + g_x_0_zzzzz_xxxxz[k];

                g_x_0_zzzzzz_xxxy[k] = -g_x_0_zzzzz_xxxy[k] * cd_z[k] + g_x_0_zzzzz_xxxyz[k];

                g_x_0_zzzzzz_xxxz[k] = -g_x_0_zzzzz_xxxz[k] * cd_z[k] + g_x_0_zzzzz_xxxzz[k];

                g_x_0_zzzzzz_xxyy[k] = -g_x_0_zzzzz_xxyy[k] * cd_z[k] + g_x_0_zzzzz_xxyyz[k];

                g_x_0_zzzzzz_xxyz[k] = -g_x_0_zzzzz_xxyz[k] * cd_z[k] + g_x_0_zzzzz_xxyzz[k];

                g_x_0_zzzzzz_xxzz[k] = -g_x_0_zzzzz_xxzz[k] * cd_z[k] + g_x_0_zzzzz_xxzzz[k];

                g_x_0_zzzzzz_xyyy[k] = -g_x_0_zzzzz_xyyy[k] * cd_z[k] + g_x_0_zzzzz_xyyyz[k];

                g_x_0_zzzzzz_xyyz[k] = -g_x_0_zzzzz_xyyz[k] * cd_z[k] + g_x_0_zzzzz_xyyzz[k];

                g_x_0_zzzzzz_xyzz[k] = -g_x_0_zzzzz_xyzz[k] * cd_z[k] + g_x_0_zzzzz_xyzzz[k];

                g_x_0_zzzzzz_xzzz[k] = -g_x_0_zzzzz_xzzz[k] * cd_z[k] + g_x_0_zzzzz_xzzzz[k];

                g_x_0_zzzzzz_yyyy[k] = -g_x_0_zzzzz_yyyy[k] * cd_z[k] + g_x_0_zzzzz_yyyyz[k];

                g_x_0_zzzzzz_yyyz[k] = -g_x_0_zzzzz_yyyz[k] * cd_z[k] + g_x_0_zzzzz_yyyzz[k];

                g_x_0_zzzzzz_yyzz[k] = -g_x_0_zzzzz_yyzz[k] * cd_z[k] + g_x_0_zzzzz_yyzzz[k];

                g_x_0_zzzzzz_yzzz[k] = -g_x_0_zzzzz_yzzz[k] * cd_z[k] + g_x_0_zzzzz_yzzzz[k];

                g_x_0_zzzzzz_zzzz[k] = -g_x_0_zzzzz_zzzz[k] * cd_z[k] + g_x_0_zzzzz_zzzzz[k];
            }
            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 0);

            auto g_y_0_xxxxxx_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 1);

            auto g_y_0_xxxxxx_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 2);

            auto g_y_0_xxxxxx_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 3);

            auto g_y_0_xxxxxx_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 4);

            auto g_y_0_xxxxxx_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 5);

            auto g_y_0_xxxxxx_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 6);

            auto g_y_0_xxxxxx_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 7);

            auto g_y_0_xxxxxx_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 8);

            auto g_y_0_xxxxxx_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 9);

            auto g_y_0_xxxxxx_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 10);

            auto g_y_0_xxxxxx_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 11);

            auto g_y_0_xxxxxx_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 12);

            auto g_y_0_xxxxxx_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 13);

            auto g_y_0_xxxxxx_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxx_xxxx, g_y_0_xxxxx_xxxxx, g_y_0_xxxxx_xxxxy, g_y_0_xxxxx_xxxxz, g_y_0_xxxxx_xxxy, g_y_0_xxxxx_xxxyy, g_y_0_xxxxx_xxxyz, g_y_0_xxxxx_xxxz, g_y_0_xxxxx_xxxzz, g_y_0_xxxxx_xxyy, g_y_0_xxxxx_xxyyy, g_y_0_xxxxx_xxyyz, g_y_0_xxxxx_xxyz, g_y_0_xxxxx_xxyzz, g_y_0_xxxxx_xxzz, g_y_0_xxxxx_xxzzz, g_y_0_xxxxx_xyyy, g_y_0_xxxxx_xyyyy, g_y_0_xxxxx_xyyyz, g_y_0_xxxxx_xyyz, g_y_0_xxxxx_xyyzz, g_y_0_xxxxx_xyzz, g_y_0_xxxxx_xyzzz, g_y_0_xxxxx_xzzz, g_y_0_xxxxx_xzzzz, g_y_0_xxxxx_yyyy, g_y_0_xxxxx_yyyz, g_y_0_xxxxx_yyzz, g_y_0_xxxxx_yzzz, g_y_0_xxxxx_zzzz, g_y_0_xxxxxx_xxxx, g_y_0_xxxxxx_xxxy, g_y_0_xxxxxx_xxxz, g_y_0_xxxxxx_xxyy, g_y_0_xxxxxx_xxyz, g_y_0_xxxxxx_xxzz, g_y_0_xxxxxx_xyyy, g_y_0_xxxxxx_xyyz, g_y_0_xxxxxx_xyzz, g_y_0_xxxxxx_xzzz, g_y_0_xxxxxx_yyyy, g_y_0_xxxxxx_yyyz, g_y_0_xxxxxx_yyzz, g_y_0_xxxxxx_yzzz, g_y_0_xxxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_xxxx[k] = -g_y_0_xxxxx_xxxx[k] * cd_x[k] + g_y_0_xxxxx_xxxxx[k];

                g_y_0_xxxxxx_xxxy[k] = -g_y_0_xxxxx_xxxy[k] * cd_x[k] + g_y_0_xxxxx_xxxxy[k];

                g_y_0_xxxxxx_xxxz[k] = -g_y_0_xxxxx_xxxz[k] * cd_x[k] + g_y_0_xxxxx_xxxxz[k];

                g_y_0_xxxxxx_xxyy[k] = -g_y_0_xxxxx_xxyy[k] * cd_x[k] + g_y_0_xxxxx_xxxyy[k];

                g_y_0_xxxxxx_xxyz[k] = -g_y_0_xxxxx_xxyz[k] * cd_x[k] + g_y_0_xxxxx_xxxyz[k];

                g_y_0_xxxxxx_xxzz[k] = -g_y_0_xxxxx_xxzz[k] * cd_x[k] + g_y_0_xxxxx_xxxzz[k];

                g_y_0_xxxxxx_xyyy[k] = -g_y_0_xxxxx_xyyy[k] * cd_x[k] + g_y_0_xxxxx_xxyyy[k];

                g_y_0_xxxxxx_xyyz[k] = -g_y_0_xxxxx_xyyz[k] * cd_x[k] + g_y_0_xxxxx_xxyyz[k];

                g_y_0_xxxxxx_xyzz[k] = -g_y_0_xxxxx_xyzz[k] * cd_x[k] + g_y_0_xxxxx_xxyzz[k];

                g_y_0_xxxxxx_xzzz[k] = -g_y_0_xxxxx_xzzz[k] * cd_x[k] + g_y_0_xxxxx_xxzzz[k];

                g_y_0_xxxxxx_yyyy[k] = -g_y_0_xxxxx_yyyy[k] * cd_x[k] + g_y_0_xxxxx_xyyyy[k];

                g_y_0_xxxxxx_yyyz[k] = -g_y_0_xxxxx_yyyz[k] * cd_x[k] + g_y_0_xxxxx_xyyyz[k];

                g_y_0_xxxxxx_yyzz[k] = -g_y_0_xxxxx_yyzz[k] * cd_x[k] + g_y_0_xxxxx_xyyzz[k];

                g_y_0_xxxxxx_yzzz[k] = -g_y_0_xxxxx_yzzz[k] * cd_x[k] + g_y_0_xxxxx_xyzzz[k];

                g_y_0_xxxxxx_zzzz[k] = -g_y_0_xxxxx_zzzz[k] * cd_x[k] + g_y_0_xxxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 15);

            auto g_y_0_xxxxxy_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 16);

            auto g_y_0_xxxxxy_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 17);

            auto g_y_0_xxxxxy_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 18);

            auto g_y_0_xxxxxy_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 19);

            auto g_y_0_xxxxxy_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 20);

            auto g_y_0_xxxxxy_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 21);

            auto g_y_0_xxxxxy_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 22);

            auto g_y_0_xxxxxy_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 23);

            auto g_y_0_xxxxxy_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 24);

            auto g_y_0_xxxxxy_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 25);

            auto g_y_0_xxxxxy_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 26);

            auto g_y_0_xxxxxy_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 27);

            auto g_y_0_xxxxxy_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 28);

            auto g_y_0_xxxxxy_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxy_xxxx, g_y_0_xxxxxy_xxxy, g_y_0_xxxxxy_xxxz, g_y_0_xxxxxy_xxyy, g_y_0_xxxxxy_xxyz, g_y_0_xxxxxy_xxzz, g_y_0_xxxxxy_xyyy, g_y_0_xxxxxy_xyyz, g_y_0_xxxxxy_xyzz, g_y_0_xxxxxy_xzzz, g_y_0_xxxxxy_yyyy, g_y_0_xxxxxy_yyyz, g_y_0_xxxxxy_yyzz, g_y_0_xxxxxy_yzzz, g_y_0_xxxxxy_zzzz, g_y_0_xxxxy_xxxx, g_y_0_xxxxy_xxxxx, g_y_0_xxxxy_xxxxy, g_y_0_xxxxy_xxxxz, g_y_0_xxxxy_xxxy, g_y_0_xxxxy_xxxyy, g_y_0_xxxxy_xxxyz, g_y_0_xxxxy_xxxz, g_y_0_xxxxy_xxxzz, g_y_0_xxxxy_xxyy, g_y_0_xxxxy_xxyyy, g_y_0_xxxxy_xxyyz, g_y_0_xxxxy_xxyz, g_y_0_xxxxy_xxyzz, g_y_0_xxxxy_xxzz, g_y_0_xxxxy_xxzzz, g_y_0_xxxxy_xyyy, g_y_0_xxxxy_xyyyy, g_y_0_xxxxy_xyyyz, g_y_0_xxxxy_xyyz, g_y_0_xxxxy_xyyzz, g_y_0_xxxxy_xyzz, g_y_0_xxxxy_xyzzz, g_y_0_xxxxy_xzzz, g_y_0_xxxxy_xzzzz, g_y_0_xxxxy_yyyy, g_y_0_xxxxy_yyyz, g_y_0_xxxxy_yyzz, g_y_0_xxxxy_yzzz, g_y_0_xxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_xxxx[k] = -g_y_0_xxxxy_xxxx[k] * cd_x[k] + g_y_0_xxxxy_xxxxx[k];

                g_y_0_xxxxxy_xxxy[k] = -g_y_0_xxxxy_xxxy[k] * cd_x[k] + g_y_0_xxxxy_xxxxy[k];

                g_y_0_xxxxxy_xxxz[k] = -g_y_0_xxxxy_xxxz[k] * cd_x[k] + g_y_0_xxxxy_xxxxz[k];

                g_y_0_xxxxxy_xxyy[k] = -g_y_0_xxxxy_xxyy[k] * cd_x[k] + g_y_0_xxxxy_xxxyy[k];

                g_y_0_xxxxxy_xxyz[k] = -g_y_0_xxxxy_xxyz[k] * cd_x[k] + g_y_0_xxxxy_xxxyz[k];

                g_y_0_xxxxxy_xxzz[k] = -g_y_0_xxxxy_xxzz[k] * cd_x[k] + g_y_0_xxxxy_xxxzz[k];

                g_y_0_xxxxxy_xyyy[k] = -g_y_0_xxxxy_xyyy[k] * cd_x[k] + g_y_0_xxxxy_xxyyy[k];

                g_y_0_xxxxxy_xyyz[k] = -g_y_0_xxxxy_xyyz[k] * cd_x[k] + g_y_0_xxxxy_xxyyz[k];

                g_y_0_xxxxxy_xyzz[k] = -g_y_0_xxxxy_xyzz[k] * cd_x[k] + g_y_0_xxxxy_xxyzz[k];

                g_y_0_xxxxxy_xzzz[k] = -g_y_0_xxxxy_xzzz[k] * cd_x[k] + g_y_0_xxxxy_xxzzz[k];

                g_y_0_xxxxxy_yyyy[k] = -g_y_0_xxxxy_yyyy[k] * cd_x[k] + g_y_0_xxxxy_xyyyy[k];

                g_y_0_xxxxxy_yyyz[k] = -g_y_0_xxxxy_yyyz[k] * cd_x[k] + g_y_0_xxxxy_xyyyz[k];

                g_y_0_xxxxxy_yyzz[k] = -g_y_0_xxxxy_yyzz[k] * cd_x[k] + g_y_0_xxxxy_xyyzz[k];

                g_y_0_xxxxxy_yzzz[k] = -g_y_0_xxxxy_yzzz[k] * cd_x[k] + g_y_0_xxxxy_xyzzz[k];

                g_y_0_xxxxxy_zzzz[k] = -g_y_0_xxxxy_zzzz[k] * cd_x[k] + g_y_0_xxxxy_xzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 30);

            auto g_y_0_xxxxxz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 31);

            auto g_y_0_xxxxxz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 32);

            auto g_y_0_xxxxxz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 33);

            auto g_y_0_xxxxxz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 34);

            auto g_y_0_xxxxxz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 35);

            auto g_y_0_xxxxxz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 36);

            auto g_y_0_xxxxxz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 37);

            auto g_y_0_xxxxxz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 38);

            auto g_y_0_xxxxxz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 39);

            auto g_y_0_xxxxxz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 40);

            auto g_y_0_xxxxxz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 41);

            auto g_y_0_xxxxxz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 42);

            auto g_y_0_xxxxxz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 43);

            auto g_y_0_xxxxxz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxz_xxxx, g_y_0_xxxxxz_xxxy, g_y_0_xxxxxz_xxxz, g_y_0_xxxxxz_xxyy, g_y_0_xxxxxz_xxyz, g_y_0_xxxxxz_xxzz, g_y_0_xxxxxz_xyyy, g_y_0_xxxxxz_xyyz, g_y_0_xxxxxz_xyzz, g_y_0_xxxxxz_xzzz, g_y_0_xxxxxz_yyyy, g_y_0_xxxxxz_yyyz, g_y_0_xxxxxz_yyzz, g_y_0_xxxxxz_yzzz, g_y_0_xxxxxz_zzzz, g_y_0_xxxxz_xxxx, g_y_0_xxxxz_xxxxx, g_y_0_xxxxz_xxxxy, g_y_0_xxxxz_xxxxz, g_y_0_xxxxz_xxxy, g_y_0_xxxxz_xxxyy, g_y_0_xxxxz_xxxyz, g_y_0_xxxxz_xxxz, g_y_0_xxxxz_xxxzz, g_y_0_xxxxz_xxyy, g_y_0_xxxxz_xxyyy, g_y_0_xxxxz_xxyyz, g_y_0_xxxxz_xxyz, g_y_0_xxxxz_xxyzz, g_y_0_xxxxz_xxzz, g_y_0_xxxxz_xxzzz, g_y_0_xxxxz_xyyy, g_y_0_xxxxz_xyyyy, g_y_0_xxxxz_xyyyz, g_y_0_xxxxz_xyyz, g_y_0_xxxxz_xyyzz, g_y_0_xxxxz_xyzz, g_y_0_xxxxz_xyzzz, g_y_0_xxxxz_xzzz, g_y_0_xxxxz_xzzzz, g_y_0_xxxxz_yyyy, g_y_0_xxxxz_yyyz, g_y_0_xxxxz_yyzz, g_y_0_xxxxz_yzzz, g_y_0_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_xxxx[k] = -g_y_0_xxxxz_xxxx[k] * cd_x[k] + g_y_0_xxxxz_xxxxx[k];

                g_y_0_xxxxxz_xxxy[k] = -g_y_0_xxxxz_xxxy[k] * cd_x[k] + g_y_0_xxxxz_xxxxy[k];

                g_y_0_xxxxxz_xxxz[k] = -g_y_0_xxxxz_xxxz[k] * cd_x[k] + g_y_0_xxxxz_xxxxz[k];

                g_y_0_xxxxxz_xxyy[k] = -g_y_0_xxxxz_xxyy[k] * cd_x[k] + g_y_0_xxxxz_xxxyy[k];

                g_y_0_xxxxxz_xxyz[k] = -g_y_0_xxxxz_xxyz[k] * cd_x[k] + g_y_0_xxxxz_xxxyz[k];

                g_y_0_xxxxxz_xxzz[k] = -g_y_0_xxxxz_xxzz[k] * cd_x[k] + g_y_0_xxxxz_xxxzz[k];

                g_y_0_xxxxxz_xyyy[k] = -g_y_0_xxxxz_xyyy[k] * cd_x[k] + g_y_0_xxxxz_xxyyy[k];

                g_y_0_xxxxxz_xyyz[k] = -g_y_0_xxxxz_xyyz[k] * cd_x[k] + g_y_0_xxxxz_xxyyz[k];

                g_y_0_xxxxxz_xyzz[k] = -g_y_0_xxxxz_xyzz[k] * cd_x[k] + g_y_0_xxxxz_xxyzz[k];

                g_y_0_xxxxxz_xzzz[k] = -g_y_0_xxxxz_xzzz[k] * cd_x[k] + g_y_0_xxxxz_xxzzz[k];

                g_y_0_xxxxxz_yyyy[k] = -g_y_0_xxxxz_yyyy[k] * cd_x[k] + g_y_0_xxxxz_xyyyy[k];

                g_y_0_xxxxxz_yyyz[k] = -g_y_0_xxxxz_yyyz[k] * cd_x[k] + g_y_0_xxxxz_xyyyz[k];

                g_y_0_xxxxxz_yyzz[k] = -g_y_0_xxxxz_yyzz[k] * cd_x[k] + g_y_0_xxxxz_xyyzz[k];

                g_y_0_xxxxxz_yzzz[k] = -g_y_0_xxxxz_yzzz[k] * cd_x[k] + g_y_0_xxxxz_xyzzz[k];

                g_y_0_xxxxxz_zzzz[k] = -g_y_0_xxxxz_zzzz[k] * cd_x[k] + g_y_0_xxxxz_xzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 45);

            auto g_y_0_xxxxyy_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 46);

            auto g_y_0_xxxxyy_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 47);

            auto g_y_0_xxxxyy_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 48);

            auto g_y_0_xxxxyy_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 49);

            auto g_y_0_xxxxyy_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 50);

            auto g_y_0_xxxxyy_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 51);

            auto g_y_0_xxxxyy_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 52);

            auto g_y_0_xxxxyy_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 53);

            auto g_y_0_xxxxyy_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 54);

            auto g_y_0_xxxxyy_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 55);

            auto g_y_0_xxxxyy_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 56);

            auto g_y_0_xxxxyy_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 57);

            auto g_y_0_xxxxyy_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 58);

            auto g_y_0_xxxxyy_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyy_xxxx, g_y_0_xxxxyy_xxxy, g_y_0_xxxxyy_xxxz, g_y_0_xxxxyy_xxyy, g_y_0_xxxxyy_xxyz, g_y_0_xxxxyy_xxzz, g_y_0_xxxxyy_xyyy, g_y_0_xxxxyy_xyyz, g_y_0_xxxxyy_xyzz, g_y_0_xxxxyy_xzzz, g_y_0_xxxxyy_yyyy, g_y_0_xxxxyy_yyyz, g_y_0_xxxxyy_yyzz, g_y_0_xxxxyy_yzzz, g_y_0_xxxxyy_zzzz, g_y_0_xxxyy_xxxx, g_y_0_xxxyy_xxxxx, g_y_0_xxxyy_xxxxy, g_y_0_xxxyy_xxxxz, g_y_0_xxxyy_xxxy, g_y_0_xxxyy_xxxyy, g_y_0_xxxyy_xxxyz, g_y_0_xxxyy_xxxz, g_y_0_xxxyy_xxxzz, g_y_0_xxxyy_xxyy, g_y_0_xxxyy_xxyyy, g_y_0_xxxyy_xxyyz, g_y_0_xxxyy_xxyz, g_y_0_xxxyy_xxyzz, g_y_0_xxxyy_xxzz, g_y_0_xxxyy_xxzzz, g_y_0_xxxyy_xyyy, g_y_0_xxxyy_xyyyy, g_y_0_xxxyy_xyyyz, g_y_0_xxxyy_xyyz, g_y_0_xxxyy_xyyzz, g_y_0_xxxyy_xyzz, g_y_0_xxxyy_xyzzz, g_y_0_xxxyy_xzzz, g_y_0_xxxyy_xzzzz, g_y_0_xxxyy_yyyy, g_y_0_xxxyy_yyyz, g_y_0_xxxyy_yyzz, g_y_0_xxxyy_yzzz, g_y_0_xxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_xxxx[k] = -g_y_0_xxxyy_xxxx[k] * cd_x[k] + g_y_0_xxxyy_xxxxx[k];

                g_y_0_xxxxyy_xxxy[k] = -g_y_0_xxxyy_xxxy[k] * cd_x[k] + g_y_0_xxxyy_xxxxy[k];

                g_y_0_xxxxyy_xxxz[k] = -g_y_0_xxxyy_xxxz[k] * cd_x[k] + g_y_0_xxxyy_xxxxz[k];

                g_y_0_xxxxyy_xxyy[k] = -g_y_0_xxxyy_xxyy[k] * cd_x[k] + g_y_0_xxxyy_xxxyy[k];

                g_y_0_xxxxyy_xxyz[k] = -g_y_0_xxxyy_xxyz[k] * cd_x[k] + g_y_0_xxxyy_xxxyz[k];

                g_y_0_xxxxyy_xxzz[k] = -g_y_0_xxxyy_xxzz[k] * cd_x[k] + g_y_0_xxxyy_xxxzz[k];

                g_y_0_xxxxyy_xyyy[k] = -g_y_0_xxxyy_xyyy[k] * cd_x[k] + g_y_0_xxxyy_xxyyy[k];

                g_y_0_xxxxyy_xyyz[k] = -g_y_0_xxxyy_xyyz[k] * cd_x[k] + g_y_0_xxxyy_xxyyz[k];

                g_y_0_xxxxyy_xyzz[k] = -g_y_0_xxxyy_xyzz[k] * cd_x[k] + g_y_0_xxxyy_xxyzz[k];

                g_y_0_xxxxyy_xzzz[k] = -g_y_0_xxxyy_xzzz[k] * cd_x[k] + g_y_0_xxxyy_xxzzz[k];

                g_y_0_xxxxyy_yyyy[k] = -g_y_0_xxxyy_yyyy[k] * cd_x[k] + g_y_0_xxxyy_xyyyy[k];

                g_y_0_xxxxyy_yyyz[k] = -g_y_0_xxxyy_yyyz[k] * cd_x[k] + g_y_0_xxxyy_xyyyz[k];

                g_y_0_xxxxyy_yyzz[k] = -g_y_0_xxxyy_yyzz[k] * cd_x[k] + g_y_0_xxxyy_xyyzz[k];

                g_y_0_xxxxyy_yzzz[k] = -g_y_0_xxxyy_yzzz[k] * cd_x[k] + g_y_0_xxxyy_xyzzz[k];

                g_y_0_xxxxyy_zzzz[k] = -g_y_0_xxxyy_zzzz[k] * cd_x[k] + g_y_0_xxxyy_xzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 60);

            auto g_y_0_xxxxyz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 61);

            auto g_y_0_xxxxyz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 62);

            auto g_y_0_xxxxyz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 63);

            auto g_y_0_xxxxyz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 64);

            auto g_y_0_xxxxyz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 65);

            auto g_y_0_xxxxyz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 66);

            auto g_y_0_xxxxyz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 67);

            auto g_y_0_xxxxyz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 68);

            auto g_y_0_xxxxyz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 69);

            auto g_y_0_xxxxyz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 70);

            auto g_y_0_xxxxyz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 71);

            auto g_y_0_xxxxyz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 72);

            auto g_y_0_xxxxyz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 73);

            auto g_y_0_xxxxyz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 74);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyz_xxxx, g_y_0_xxxxyz_xxxy, g_y_0_xxxxyz_xxxz, g_y_0_xxxxyz_xxyy, g_y_0_xxxxyz_xxyz, g_y_0_xxxxyz_xxzz, g_y_0_xxxxyz_xyyy, g_y_0_xxxxyz_xyyz, g_y_0_xxxxyz_xyzz, g_y_0_xxxxyz_xzzz, g_y_0_xxxxyz_yyyy, g_y_0_xxxxyz_yyyz, g_y_0_xxxxyz_yyzz, g_y_0_xxxxyz_yzzz, g_y_0_xxxxyz_zzzz, g_y_0_xxxyz_xxxx, g_y_0_xxxyz_xxxxx, g_y_0_xxxyz_xxxxy, g_y_0_xxxyz_xxxxz, g_y_0_xxxyz_xxxy, g_y_0_xxxyz_xxxyy, g_y_0_xxxyz_xxxyz, g_y_0_xxxyz_xxxz, g_y_0_xxxyz_xxxzz, g_y_0_xxxyz_xxyy, g_y_0_xxxyz_xxyyy, g_y_0_xxxyz_xxyyz, g_y_0_xxxyz_xxyz, g_y_0_xxxyz_xxyzz, g_y_0_xxxyz_xxzz, g_y_0_xxxyz_xxzzz, g_y_0_xxxyz_xyyy, g_y_0_xxxyz_xyyyy, g_y_0_xxxyz_xyyyz, g_y_0_xxxyz_xyyz, g_y_0_xxxyz_xyyzz, g_y_0_xxxyz_xyzz, g_y_0_xxxyz_xyzzz, g_y_0_xxxyz_xzzz, g_y_0_xxxyz_xzzzz, g_y_0_xxxyz_yyyy, g_y_0_xxxyz_yyyz, g_y_0_xxxyz_yyzz, g_y_0_xxxyz_yzzz, g_y_0_xxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_xxxx[k] = -g_y_0_xxxyz_xxxx[k] * cd_x[k] + g_y_0_xxxyz_xxxxx[k];

                g_y_0_xxxxyz_xxxy[k] = -g_y_0_xxxyz_xxxy[k] * cd_x[k] + g_y_0_xxxyz_xxxxy[k];

                g_y_0_xxxxyz_xxxz[k] = -g_y_0_xxxyz_xxxz[k] * cd_x[k] + g_y_0_xxxyz_xxxxz[k];

                g_y_0_xxxxyz_xxyy[k] = -g_y_0_xxxyz_xxyy[k] * cd_x[k] + g_y_0_xxxyz_xxxyy[k];

                g_y_0_xxxxyz_xxyz[k] = -g_y_0_xxxyz_xxyz[k] * cd_x[k] + g_y_0_xxxyz_xxxyz[k];

                g_y_0_xxxxyz_xxzz[k] = -g_y_0_xxxyz_xxzz[k] * cd_x[k] + g_y_0_xxxyz_xxxzz[k];

                g_y_0_xxxxyz_xyyy[k] = -g_y_0_xxxyz_xyyy[k] * cd_x[k] + g_y_0_xxxyz_xxyyy[k];

                g_y_0_xxxxyz_xyyz[k] = -g_y_0_xxxyz_xyyz[k] * cd_x[k] + g_y_0_xxxyz_xxyyz[k];

                g_y_0_xxxxyz_xyzz[k] = -g_y_0_xxxyz_xyzz[k] * cd_x[k] + g_y_0_xxxyz_xxyzz[k];

                g_y_0_xxxxyz_xzzz[k] = -g_y_0_xxxyz_xzzz[k] * cd_x[k] + g_y_0_xxxyz_xxzzz[k];

                g_y_0_xxxxyz_yyyy[k] = -g_y_0_xxxyz_yyyy[k] * cd_x[k] + g_y_0_xxxyz_xyyyy[k];

                g_y_0_xxxxyz_yyyz[k] = -g_y_0_xxxyz_yyyz[k] * cd_x[k] + g_y_0_xxxyz_xyyyz[k];

                g_y_0_xxxxyz_yyzz[k] = -g_y_0_xxxyz_yyzz[k] * cd_x[k] + g_y_0_xxxyz_xyyzz[k];

                g_y_0_xxxxyz_yzzz[k] = -g_y_0_xxxyz_yzzz[k] * cd_x[k] + g_y_0_xxxyz_xyzzz[k];

                g_y_0_xxxxyz_zzzz[k] = -g_y_0_xxxyz_zzzz[k] * cd_x[k] + g_y_0_xxxyz_xzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 75);

            auto g_y_0_xxxxzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 76);

            auto g_y_0_xxxxzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 77);

            auto g_y_0_xxxxzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 78);

            auto g_y_0_xxxxzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 79);

            auto g_y_0_xxxxzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 80);

            auto g_y_0_xxxxzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 81);

            auto g_y_0_xxxxzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 82);

            auto g_y_0_xxxxzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 83);

            auto g_y_0_xxxxzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 84);

            auto g_y_0_xxxxzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 85);

            auto g_y_0_xxxxzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 86);

            auto g_y_0_xxxxzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 87);

            auto g_y_0_xxxxzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 88);

            auto g_y_0_xxxxzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxzz_xxxx, g_y_0_xxxxzz_xxxy, g_y_0_xxxxzz_xxxz, g_y_0_xxxxzz_xxyy, g_y_0_xxxxzz_xxyz, g_y_0_xxxxzz_xxzz, g_y_0_xxxxzz_xyyy, g_y_0_xxxxzz_xyyz, g_y_0_xxxxzz_xyzz, g_y_0_xxxxzz_xzzz, g_y_0_xxxxzz_yyyy, g_y_0_xxxxzz_yyyz, g_y_0_xxxxzz_yyzz, g_y_0_xxxxzz_yzzz, g_y_0_xxxxzz_zzzz, g_y_0_xxxzz_xxxx, g_y_0_xxxzz_xxxxx, g_y_0_xxxzz_xxxxy, g_y_0_xxxzz_xxxxz, g_y_0_xxxzz_xxxy, g_y_0_xxxzz_xxxyy, g_y_0_xxxzz_xxxyz, g_y_0_xxxzz_xxxz, g_y_0_xxxzz_xxxzz, g_y_0_xxxzz_xxyy, g_y_0_xxxzz_xxyyy, g_y_0_xxxzz_xxyyz, g_y_0_xxxzz_xxyz, g_y_0_xxxzz_xxyzz, g_y_0_xxxzz_xxzz, g_y_0_xxxzz_xxzzz, g_y_0_xxxzz_xyyy, g_y_0_xxxzz_xyyyy, g_y_0_xxxzz_xyyyz, g_y_0_xxxzz_xyyz, g_y_0_xxxzz_xyyzz, g_y_0_xxxzz_xyzz, g_y_0_xxxzz_xyzzz, g_y_0_xxxzz_xzzz, g_y_0_xxxzz_xzzzz, g_y_0_xxxzz_yyyy, g_y_0_xxxzz_yyyz, g_y_0_xxxzz_yyzz, g_y_0_xxxzz_yzzz, g_y_0_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_xxxx[k] = -g_y_0_xxxzz_xxxx[k] * cd_x[k] + g_y_0_xxxzz_xxxxx[k];

                g_y_0_xxxxzz_xxxy[k] = -g_y_0_xxxzz_xxxy[k] * cd_x[k] + g_y_0_xxxzz_xxxxy[k];

                g_y_0_xxxxzz_xxxz[k] = -g_y_0_xxxzz_xxxz[k] * cd_x[k] + g_y_0_xxxzz_xxxxz[k];

                g_y_0_xxxxzz_xxyy[k] = -g_y_0_xxxzz_xxyy[k] * cd_x[k] + g_y_0_xxxzz_xxxyy[k];

                g_y_0_xxxxzz_xxyz[k] = -g_y_0_xxxzz_xxyz[k] * cd_x[k] + g_y_0_xxxzz_xxxyz[k];

                g_y_0_xxxxzz_xxzz[k] = -g_y_0_xxxzz_xxzz[k] * cd_x[k] + g_y_0_xxxzz_xxxzz[k];

                g_y_0_xxxxzz_xyyy[k] = -g_y_0_xxxzz_xyyy[k] * cd_x[k] + g_y_0_xxxzz_xxyyy[k];

                g_y_0_xxxxzz_xyyz[k] = -g_y_0_xxxzz_xyyz[k] * cd_x[k] + g_y_0_xxxzz_xxyyz[k];

                g_y_0_xxxxzz_xyzz[k] = -g_y_0_xxxzz_xyzz[k] * cd_x[k] + g_y_0_xxxzz_xxyzz[k];

                g_y_0_xxxxzz_xzzz[k] = -g_y_0_xxxzz_xzzz[k] * cd_x[k] + g_y_0_xxxzz_xxzzz[k];

                g_y_0_xxxxzz_yyyy[k] = -g_y_0_xxxzz_yyyy[k] * cd_x[k] + g_y_0_xxxzz_xyyyy[k];

                g_y_0_xxxxzz_yyyz[k] = -g_y_0_xxxzz_yyyz[k] * cd_x[k] + g_y_0_xxxzz_xyyyz[k];

                g_y_0_xxxxzz_yyzz[k] = -g_y_0_xxxzz_yyzz[k] * cd_x[k] + g_y_0_xxxzz_xyyzz[k];

                g_y_0_xxxxzz_yzzz[k] = -g_y_0_xxxzz_yzzz[k] * cd_x[k] + g_y_0_xxxzz_xyzzz[k];

                g_y_0_xxxxzz_zzzz[k] = -g_y_0_xxxzz_zzzz[k] * cd_x[k] + g_y_0_xxxzz_xzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 90);

            auto g_y_0_xxxyyy_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 91);

            auto g_y_0_xxxyyy_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 92);

            auto g_y_0_xxxyyy_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 93);

            auto g_y_0_xxxyyy_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 94);

            auto g_y_0_xxxyyy_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 95);

            auto g_y_0_xxxyyy_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 96);

            auto g_y_0_xxxyyy_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 97);

            auto g_y_0_xxxyyy_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 98);

            auto g_y_0_xxxyyy_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 99);

            auto g_y_0_xxxyyy_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 100);

            auto g_y_0_xxxyyy_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 101);

            auto g_y_0_xxxyyy_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 102);

            auto g_y_0_xxxyyy_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 103);

            auto g_y_0_xxxyyy_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyy_xxxx, g_y_0_xxxyyy_xxxy, g_y_0_xxxyyy_xxxz, g_y_0_xxxyyy_xxyy, g_y_0_xxxyyy_xxyz, g_y_0_xxxyyy_xxzz, g_y_0_xxxyyy_xyyy, g_y_0_xxxyyy_xyyz, g_y_0_xxxyyy_xyzz, g_y_0_xxxyyy_xzzz, g_y_0_xxxyyy_yyyy, g_y_0_xxxyyy_yyyz, g_y_0_xxxyyy_yyzz, g_y_0_xxxyyy_yzzz, g_y_0_xxxyyy_zzzz, g_y_0_xxyyy_xxxx, g_y_0_xxyyy_xxxxx, g_y_0_xxyyy_xxxxy, g_y_0_xxyyy_xxxxz, g_y_0_xxyyy_xxxy, g_y_0_xxyyy_xxxyy, g_y_0_xxyyy_xxxyz, g_y_0_xxyyy_xxxz, g_y_0_xxyyy_xxxzz, g_y_0_xxyyy_xxyy, g_y_0_xxyyy_xxyyy, g_y_0_xxyyy_xxyyz, g_y_0_xxyyy_xxyz, g_y_0_xxyyy_xxyzz, g_y_0_xxyyy_xxzz, g_y_0_xxyyy_xxzzz, g_y_0_xxyyy_xyyy, g_y_0_xxyyy_xyyyy, g_y_0_xxyyy_xyyyz, g_y_0_xxyyy_xyyz, g_y_0_xxyyy_xyyzz, g_y_0_xxyyy_xyzz, g_y_0_xxyyy_xyzzz, g_y_0_xxyyy_xzzz, g_y_0_xxyyy_xzzzz, g_y_0_xxyyy_yyyy, g_y_0_xxyyy_yyyz, g_y_0_xxyyy_yyzz, g_y_0_xxyyy_yzzz, g_y_0_xxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_xxxx[k] = -g_y_0_xxyyy_xxxx[k] * cd_x[k] + g_y_0_xxyyy_xxxxx[k];

                g_y_0_xxxyyy_xxxy[k] = -g_y_0_xxyyy_xxxy[k] * cd_x[k] + g_y_0_xxyyy_xxxxy[k];

                g_y_0_xxxyyy_xxxz[k] = -g_y_0_xxyyy_xxxz[k] * cd_x[k] + g_y_0_xxyyy_xxxxz[k];

                g_y_0_xxxyyy_xxyy[k] = -g_y_0_xxyyy_xxyy[k] * cd_x[k] + g_y_0_xxyyy_xxxyy[k];

                g_y_0_xxxyyy_xxyz[k] = -g_y_0_xxyyy_xxyz[k] * cd_x[k] + g_y_0_xxyyy_xxxyz[k];

                g_y_0_xxxyyy_xxzz[k] = -g_y_0_xxyyy_xxzz[k] * cd_x[k] + g_y_0_xxyyy_xxxzz[k];

                g_y_0_xxxyyy_xyyy[k] = -g_y_0_xxyyy_xyyy[k] * cd_x[k] + g_y_0_xxyyy_xxyyy[k];

                g_y_0_xxxyyy_xyyz[k] = -g_y_0_xxyyy_xyyz[k] * cd_x[k] + g_y_0_xxyyy_xxyyz[k];

                g_y_0_xxxyyy_xyzz[k] = -g_y_0_xxyyy_xyzz[k] * cd_x[k] + g_y_0_xxyyy_xxyzz[k];

                g_y_0_xxxyyy_xzzz[k] = -g_y_0_xxyyy_xzzz[k] * cd_x[k] + g_y_0_xxyyy_xxzzz[k];

                g_y_0_xxxyyy_yyyy[k] = -g_y_0_xxyyy_yyyy[k] * cd_x[k] + g_y_0_xxyyy_xyyyy[k];

                g_y_0_xxxyyy_yyyz[k] = -g_y_0_xxyyy_yyyz[k] * cd_x[k] + g_y_0_xxyyy_xyyyz[k];

                g_y_0_xxxyyy_yyzz[k] = -g_y_0_xxyyy_yyzz[k] * cd_x[k] + g_y_0_xxyyy_xyyzz[k];

                g_y_0_xxxyyy_yzzz[k] = -g_y_0_xxyyy_yzzz[k] * cd_x[k] + g_y_0_xxyyy_xyzzz[k];

                g_y_0_xxxyyy_zzzz[k] = -g_y_0_xxyyy_zzzz[k] * cd_x[k] + g_y_0_xxyyy_xzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 105);

            auto g_y_0_xxxyyz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 106);

            auto g_y_0_xxxyyz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 107);

            auto g_y_0_xxxyyz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 108);

            auto g_y_0_xxxyyz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 109);

            auto g_y_0_xxxyyz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 110);

            auto g_y_0_xxxyyz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 111);

            auto g_y_0_xxxyyz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 112);

            auto g_y_0_xxxyyz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 113);

            auto g_y_0_xxxyyz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 114);

            auto g_y_0_xxxyyz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 115);

            auto g_y_0_xxxyyz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 116);

            auto g_y_0_xxxyyz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 117);

            auto g_y_0_xxxyyz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 118);

            auto g_y_0_xxxyyz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyz_xxxx, g_y_0_xxxyyz_xxxy, g_y_0_xxxyyz_xxxz, g_y_0_xxxyyz_xxyy, g_y_0_xxxyyz_xxyz, g_y_0_xxxyyz_xxzz, g_y_0_xxxyyz_xyyy, g_y_0_xxxyyz_xyyz, g_y_0_xxxyyz_xyzz, g_y_0_xxxyyz_xzzz, g_y_0_xxxyyz_yyyy, g_y_0_xxxyyz_yyyz, g_y_0_xxxyyz_yyzz, g_y_0_xxxyyz_yzzz, g_y_0_xxxyyz_zzzz, g_y_0_xxyyz_xxxx, g_y_0_xxyyz_xxxxx, g_y_0_xxyyz_xxxxy, g_y_0_xxyyz_xxxxz, g_y_0_xxyyz_xxxy, g_y_0_xxyyz_xxxyy, g_y_0_xxyyz_xxxyz, g_y_0_xxyyz_xxxz, g_y_0_xxyyz_xxxzz, g_y_0_xxyyz_xxyy, g_y_0_xxyyz_xxyyy, g_y_0_xxyyz_xxyyz, g_y_0_xxyyz_xxyz, g_y_0_xxyyz_xxyzz, g_y_0_xxyyz_xxzz, g_y_0_xxyyz_xxzzz, g_y_0_xxyyz_xyyy, g_y_0_xxyyz_xyyyy, g_y_0_xxyyz_xyyyz, g_y_0_xxyyz_xyyz, g_y_0_xxyyz_xyyzz, g_y_0_xxyyz_xyzz, g_y_0_xxyyz_xyzzz, g_y_0_xxyyz_xzzz, g_y_0_xxyyz_xzzzz, g_y_0_xxyyz_yyyy, g_y_0_xxyyz_yyyz, g_y_0_xxyyz_yyzz, g_y_0_xxyyz_yzzz, g_y_0_xxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_xxxx[k] = -g_y_0_xxyyz_xxxx[k] * cd_x[k] + g_y_0_xxyyz_xxxxx[k];

                g_y_0_xxxyyz_xxxy[k] = -g_y_0_xxyyz_xxxy[k] * cd_x[k] + g_y_0_xxyyz_xxxxy[k];

                g_y_0_xxxyyz_xxxz[k] = -g_y_0_xxyyz_xxxz[k] * cd_x[k] + g_y_0_xxyyz_xxxxz[k];

                g_y_0_xxxyyz_xxyy[k] = -g_y_0_xxyyz_xxyy[k] * cd_x[k] + g_y_0_xxyyz_xxxyy[k];

                g_y_0_xxxyyz_xxyz[k] = -g_y_0_xxyyz_xxyz[k] * cd_x[k] + g_y_0_xxyyz_xxxyz[k];

                g_y_0_xxxyyz_xxzz[k] = -g_y_0_xxyyz_xxzz[k] * cd_x[k] + g_y_0_xxyyz_xxxzz[k];

                g_y_0_xxxyyz_xyyy[k] = -g_y_0_xxyyz_xyyy[k] * cd_x[k] + g_y_0_xxyyz_xxyyy[k];

                g_y_0_xxxyyz_xyyz[k] = -g_y_0_xxyyz_xyyz[k] * cd_x[k] + g_y_0_xxyyz_xxyyz[k];

                g_y_0_xxxyyz_xyzz[k] = -g_y_0_xxyyz_xyzz[k] * cd_x[k] + g_y_0_xxyyz_xxyzz[k];

                g_y_0_xxxyyz_xzzz[k] = -g_y_0_xxyyz_xzzz[k] * cd_x[k] + g_y_0_xxyyz_xxzzz[k];

                g_y_0_xxxyyz_yyyy[k] = -g_y_0_xxyyz_yyyy[k] * cd_x[k] + g_y_0_xxyyz_xyyyy[k];

                g_y_0_xxxyyz_yyyz[k] = -g_y_0_xxyyz_yyyz[k] * cd_x[k] + g_y_0_xxyyz_xyyyz[k];

                g_y_0_xxxyyz_yyzz[k] = -g_y_0_xxyyz_yyzz[k] * cd_x[k] + g_y_0_xxyyz_xyyzz[k];

                g_y_0_xxxyyz_yzzz[k] = -g_y_0_xxyyz_yzzz[k] * cd_x[k] + g_y_0_xxyyz_xyzzz[k];

                g_y_0_xxxyyz_zzzz[k] = -g_y_0_xxyyz_zzzz[k] * cd_x[k] + g_y_0_xxyyz_xzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 120);

            auto g_y_0_xxxyzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 121);

            auto g_y_0_xxxyzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 122);

            auto g_y_0_xxxyzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 123);

            auto g_y_0_xxxyzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 124);

            auto g_y_0_xxxyzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 125);

            auto g_y_0_xxxyzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 126);

            auto g_y_0_xxxyzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 127);

            auto g_y_0_xxxyzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 128);

            auto g_y_0_xxxyzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 129);

            auto g_y_0_xxxyzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 130);

            auto g_y_0_xxxyzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 131);

            auto g_y_0_xxxyzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 132);

            auto g_y_0_xxxyzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 133);

            auto g_y_0_xxxyzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 134);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyzz_xxxx, g_y_0_xxxyzz_xxxy, g_y_0_xxxyzz_xxxz, g_y_0_xxxyzz_xxyy, g_y_0_xxxyzz_xxyz, g_y_0_xxxyzz_xxzz, g_y_0_xxxyzz_xyyy, g_y_0_xxxyzz_xyyz, g_y_0_xxxyzz_xyzz, g_y_0_xxxyzz_xzzz, g_y_0_xxxyzz_yyyy, g_y_0_xxxyzz_yyyz, g_y_0_xxxyzz_yyzz, g_y_0_xxxyzz_yzzz, g_y_0_xxxyzz_zzzz, g_y_0_xxyzz_xxxx, g_y_0_xxyzz_xxxxx, g_y_0_xxyzz_xxxxy, g_y_0_xxyzz_xxxxz, g_y_0_xxyzz_xxxy, g_y_0_xxyzz_xxxyy, g_y_0_xxyzz_xxxyz, g_y_0_xxyzz_xxxz, g_y_0_xxyzz_xxxzz, g_y_0_xxyzz_xxyy, g_y_0_xxyzz_xxyyy, g_y_0_xxyzz_xxyyz, g_y_0_xxyzz_xxyz, g_y_0_xxyzz_xxyzz, g_y_0_xxyzz_xxzz, g_y_0_xxyzz_xxzzz, g_y_0_xxyzz_xyyy, g_y_0_xxyzz_xyyyy, g_y_0_xxyzz_xyyyz, g_y_0_xxyzz_xyyz, g_y_0_xxyzz_xyyzz, g_y_0_xxyzz_xyzz, g_y_0_xxyzz_xyzzz, g_y_0_xxyzz_xzzz, g_y_0_xxyzz_xzzzz, g_y_0_xxyzz_yyyy, g_y_0_xxyzz_yyyz, g_y_0_xxyzz_yyzz, g_y_0_xxyzz_yzzz, g_y_0_xxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_xxxx[k] = -g_y_0_xxyzz_xxxx[k] * cd_x[k] + g_y_0_xxyzz_xxxxx[k];

                g_y_0_xxxyzz_xxxy[k] = -g_y_0_xxyzz_xxxy[k] * cd_x[k] + g_y_0_xxyzz_xxxxy[k];

                g_y_0_xxxyzz_xxxz[k] = -g_y_0_xxyzz_xxxz[k] * cd_x[k] + g_y_0_xxyzz_xxxxz[k];

                g_y_0_xxxyzz_xxyy[k] = -g_y_0_xxyzz_xxyy[k] * cd_x[k] + g_y_0_xxyzz_xxxyy[k];

                g_y_0_xxxyzz_xxyz[k] = -g_y_0_xxyzz_xxyz[k] * cd_x[k] + g_y_0_xxyzz_xxxyz[k];

                g_y_0_xxxyzz_xxzz[k] = -g_y_0_xxyzz_xxzz[k] * cd_x[k] + g_y_0_xxyzz_xxxzz[k];

                g_y_0_xxxyzz_xyyy[k] = -g_y_0_xxyzz_xyyy[k] * cd_x[k] + g_y_0_xxyzz_xxyyy[k];

                g_y_0_xxxyzz_xyyz[k] = -g_y_0_xxyzz_xyyz[k] * cd_x[k] + g_y_0_xxyzz_xxyyz[k];

                g_y_0_xxxyzz_xyzz[k] = -g_y_0_xxyzz_xyzz[k] * cd_x[k] + g_y_0_xxyzz_xxyzz[k];

                g_y_0_xxxyzz_xzzz[k] = -g_y_0_xxyzz_xzzz[k] * cd_x[k] + g_y_0_xxyzz_xxzzz[k];

                g_y_0_xxxyzz_yyyy[k] = -g_y_0_xxyzz_yyyy[k] * cd_x[k] + g_y_0_xxyzz_xyyyy[k];

                g_y_0_xxxyzz_yyyz[k] = -g_y_0_xxyzz_yyyz[k] * cd_x[k] + g_y_0_xxyzz_xyyyz[k];

                g_y_0_xxxyzz_yyzz[k] = -g_y_0_xxyzz_yyzz[k] * cd_x[k] + g_y_0_xxyzz_xyyzz[k];

                g_y_0_xxxyzz_yzzz[k] = -g_y_0_xxyzz_yzzz[k] * cd_x[k] + g_y_0_xxyzz_xyzzz[k];

                g_y_0_xxxyzz_zzzz[k] = -g_y_0_xxyzz_zzzz[k] * cd_x[k] + g_y_0_xxyzz_xzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 135);

            auto g_y_0_xxxzzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 136);

            auto g_y_0_xxxzzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 137);

            auto g_y_0_xxxzzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 138);

            auto g_y_0_xxxzzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 139);

            auto g_y_0_xxxzzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 140);

            auto g_y_0_xxxzzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 141);

            auto g_y_0_xxxzzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 142);

            auto g_y_0_xxxzzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 143);

            auto g_y_0_xxxzzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 144);

            auto g_y_0_xxxzzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 145);

            auto g_y_0_xxxzzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 146);

            auto g_y_0_xxxzzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 147);

            auto g_y_0_xxxzzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 148);

            auto g_y_0_xxxzzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzzz_xxxx, g_y_0_xxxzzz_xxxy, g_y_0_xxxzzz_xxxz, g_y_0_xxxzzz_xxyy, g_y_0_xxxzzz_xxyz, g_y_0_xxxzzz_xxzz, g_y_0_xxxzzz_xyyy, g_y_0_xxxzzz_xyyz, g_y_0_xxxzzz_xyzz, g_y_0_xxxzzz_xzzz, g_y_0_xxxzzz_yyyy, g_y_0_xxxzzz_yyyz, g_y_0_xxxzzz_yyzz, g_y_0_xxxzzz_yzzz, g_y_0_xxxzzz_zzzz, g_y_0_xxzzz_xxxx, g_y_0_xxzzz_xxxxx, g_y_0_xxzzz_xxxxy, g_y_0_xxzzz_xxxxz, g_y_0_xxzzz_xxxy, g_y_0_xxzzz_xxxyy, g_y_0_xxzzz_xxxyz, g_y_0_xxzzz_xxxz, g_y_0_xxzzz_xxxzz, g_y_0_xxzzz_xxyy, g_y_0_xxzzz_xxyyy, g_y_0_xxzzz_xxyyz, g_y_0_xxzzz_xxyz, g_y_0_xxzzz_xxyzz, g_y_0_xxzzz_xxzz, g_y_0_xxzzz_xxzzz, g_y_0_xxzzz_xyyy, g_y_0_xxzzz_xyyyy, g_y_0_xxzzz_xyyyz, g_y_0_xxzzz_xyyz, g_y_0_xxzzz_xyyzz, g_y_0_xxzzz_xyzz, g_y_0_xxzzz_xyzzz, g_y_0_xxzzz_xzzz, g_y_0_xxzzz_xzzzz, g_y_0_xxzzz_yyyy, g_y_0_xxzzz_yyyz, g_y_0_xxzzz_yyzz, g_y_0_xxzzz_yzzz, g_y_0_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_xxxx[k] = -g_y_0_xxzzz_xxxx[k] * cd_x[k] + g_y_0_xxzzz_xxxxx[k];

                g_y_0_xxxzzz_xxxy[k] = -g_y_0_xxzzz_xxxy[k] * cd_x[k] + g_y_0_xxzzz_xxxxy[k];

                g_y_0_xxxzzz_xxxz[k] = -g_y_0_xxzzz_xxxz[k] * cd_x[k] + g_y_0_xxzzz_xxxxz[k];

                g_y_0_xxxzzz_xxyy[k] = -g_y_0_xxzzz_xxyy[k] * cd_x[k] + g_y_0_xxzzz_xxxyy[k];

                g_y_0_xxxzzz_xxyz[k] = -g_y_0_xxzzz_xxyz[k] * cd_x[k] + g_y_0_xxzzz_xxxyz[k];

                g_y_0_xxxzzz_xxzz[k] = -g_y_0_xxzzz_xxzz[k] * cd_x[k] + g_y_0_xxzzz_xxxzz[k];

                g_y_0_xxxzzz_xyyy[k] = -g_y_0_xxzzz_xyyy[k] * cd_x[k] + g_y_0_xxzzz_xxyyy[k];

                g_y_0_xxxzzz_xyyz[k] = -g_y_0_xxzzz_xyyz[k] * cd_x[k] + g_y_0_xxzzz_xxyyz[k];

                g_y_0_xxxzzz_xyzz[k] = -g_y_0_xxzzz_xyzz[k] * cd_x[k] + g_y_0_xxzzz_xxyzz[k];

                g_y_0_xxxzzz_xzzz[k] = -g_y_0_xxzzz_xzzz[k] * cd_x[k] + g_y_0_xxzzz_xxzzz[k];

                g_y_0_xxxzzz_yyyy[k] = -g_y_0_xxzzz_yyyy[k] * cd_x[k] + g_y_0_xxzzz_xyyyy[k];

                g_y_0_xxxzzz_yyyz[k] = -g_y_0_xxzzz_yyyz[k] * cd_x[k] + g_y_0_xxzzz_xyyyz[k];

                g_y_0_xxxzzz_yyzz[k] = -g_y_0_xxzzz_yyzz[k] * cd_x[k] + g_y_0_xxzzz_xyyzz[k];

                g_y_0_xxxzzz_yzzz[k] = -g_y_0_xxzzz_yzzz[k] * cd_x[k] + g_y_0_xxzzz_xyzzz[k];

                g_y_0_xxxzzz_zzzz[k] = -g_y_0_xxzzz_zzzz[k] * cd_x[k] + g_y_0_xxzzz_xzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 150);

            auto g_y_0_xxyyyy_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 151);

            auto g_y_0_xxyyyy_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 152);

            auto g_y_0_xxyyyy_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 153);

            auto g_y_0_xxyyyy_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 154);

            auto g_y_0_xxyyyy_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 155);

            auto g_y_0_xxyyyy_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 156);

            auto g_y_0_xxyyyy_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 157);

            auto g_y_0_xxyyyy_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 158);

            auto g_y_0_xxyyyy_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 159);

            auto g_y_0_xxyyyy_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 160);

            auto g_y_0_xxyyyy_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 161);

            auto g_y_0_xxyyyy_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 162);

            auto g_y_0_xxyyyy_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 163);

            auto g_y_0_xxyyyy_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 164);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyy_xxxx, g_y_0_xxyyyy_xxxy, g_y_0_xxyyyy_xxxz, g_y_0_xxyyyy_xxyy, g_y_0_xxyyyy_xxyz, g_y_0_xxyyyy_xxzz, g_y_0_xxyyyy_xyyy, g_y_0_xxyyyy_xyyz, g_y_0_xxyyyy_xyzz, g_y_0_xxyyyy_xzzz, g_y_0_xxyyyy_yyyy, g_y_0_xxyyyy_yyyz, g_y_0_xxyyyy_yyzz, g_y_0_xxyyyy_yzzz, g_y_0_xxyyyy_zzzz, g_y_0_xyyyy_xxxx, g_y_0_xyyyy_xxxxx, g_y_0_xyyyy_xxxxy, g_y_0_xyyyy_xxxxz, g_y_0_xyyyy_xxxy, g_y_0_xyyyy_xxxyy, g_y_0_xyyyy_xxxyz, g_y_0_xyyyy_xxxz, g_y_0_xyyyy_xxxzz, g_y_0_xyyyy_xxyy, g_y_0_xyyyy_xxyyy, g_y_0_xyyyy_xxyyz, g_y_0_xyyyy_xxyz, g_y_0_xyyyy_xxyzz, g_y_0_xyyyy_xxzz, g_y_0_xyyyy_xxzzz, g_y_0_xyyyy_xyyy, g_y_0_xyyyy_xyyyy, g_y_0_xyyyy_xyyyz, g_y_0_xyyyy_xyyz, g_y_0_xyyyy_xyyzz, g_y_0_xyyyy_xyzz, g_y_0_xyyyy_xyzzz, g_y_0_xyyyy_xzzz, g_y_0_xyyyy_xzzzz, g_y_0_xyyyy_yyyy, g_y_0_xyyyy_yyyz, g_y_0_xyyyy_yyzz, g_y_0_xyyyy_yzzz, g_y_0_xyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_xxxx[k] = -g_y_0_xyyyy_xxxx[k] * cd_x[k] + g_y_0_xyyyy_xxxxx[k];

                g_y_0_xxyyyy_xxxy[k] = -g_y_0_xyyyy_xxxy[k] * cd_x[k] + g_y_0_xyyyy_xxxxy[k];

                g_y_0_xxyyyy_xxxz[k] = -g_y_0_xyyyy_xxxz[k] * cd_x[k] + g_y_0_xyyyy_xxxxz[k];

                g_y_0_xxyyyy_xxyy[k] = -g_y_0_xyyyy_xxyy[k] * cd_x[k] + g_y_0_xyyyy_xxxyy[k];

                g_y_0_xxyyyy_xxyz[k] = -g_y_0_xyyyy_xxyz[k] * cd_x[k] + g_y_0_xyyyy_xxxyz[k];

                g_y_0_xxyyyy_xxzz[k] = -g_y_0_xyyyy_xxzz[k] * cd_x[k] + g_y_0_xyyyy_xxxzz[k];

                g_y_0_xxyyyy_xyyy[k] = -g_y_0_xyyyy_xyyy[k] * cd_x[k] + g_y_0_xyyyy_xxyyy[k];

                g_y_0_xxyyyy_xyyz[k] = -g_y_0_xyyyy_xyyz[k] * cd_x[k] + g_y_0_xyyyy_xxyyz[k];

                g_y_0_xxyyyy_xyzz[k] = -g_y_0_xyyyy_xyzz[k] * cd_x[k] + g_y_0_xyyyy_xxyzz[k];

                g_y_0_xxyyyy_xzzz[k] = -g_y_0_xyyyy_xzzz[k] * cd_x[k] + g_y_0_xyyyy_xxzzz[k];

                g_y_0_xxyyyy_yyyy[k] = -g_y_0_xyyyy_yyyy[k] * cd_x[k] + g_y_0_xyyyy_xyyyy[k];

                g_y_0_xxyyyy_yyyz[k] = -g_y_0_xyyyy_yyyz[k] * cd_x[k] + g_y_0_xyyyy_xyyyz[k];

                g_y_0_xxyyyy_yyzz[k] = -g_y_0_xyyyy_yyzz[k] * cd_x[k] + g_y_0_xyyyy_xyyzz[k];

                g_y_0_xxyyyy_yzzz[k] = -g_y_0_xyyyy_yzzz[k] * cd_x[k] + g_y_0_xyyyy_xyzzz[k];

                g_y_0_xxyyyy_zzzz[k] = -g_y_0_xyyyy_zzzz[k] * cd_x[k] + g_y_0_xyyyy_xzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 165);

            auto g_y_0_xxyyyz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 166);

            auto g_y_0_xxyyyz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 167);

            auto g_y_0_xxyyyz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 168);

            auto g_y_0_xxyyyz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 169);

            auto g_y_0_xxyyyz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 170);

            auto g_y_0_xxyyyz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 171);

            auto g_y_0_xxyyyz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 172);

            auto g_y_0_xxyyyz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 173);

            auto g_y_0_xxyyyz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 174);

            auto g_y_0_xxyyyz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 175);

            auto g_y_0_xxyyyz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 176);

            auto g_y_0_xxyyyz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 177);

            auto g_y_0_xxyyyz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 178);

            auto g_y_0_xxyyyz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyz_xxxx, g_y_0_xxyyyz_xxxy, g_y_0_xxyyyz_xxxz, g_y_0_xxyyyz_xxyy, g_y_0_xxyyyz_xxyz, g_y_0_xxyyyz_xxzz, g_y_0_xxyyyz_xyyy, g_y_0_xxyyyz_xyyz, g_y_0_xxyyyz_xyzz, g_y_0_xxyyyz_xzzz, g_y_0_xxyyyz_yyyy, g_y_0_xxyyyz_yyyz, g_y_0_xxyyyz_yyzz, g_y_0_xxyyyz_yzzz, g_y_0_xxyyyz_zzzz, g_y_0_xyyyz_xxxx, g_y_0_xyyyz_xxxxx, g_y_0_xyyyz_xxxxy, g_y_0_xyyyz_xxxxz, g_y_0_xyyyz_xxxy, g_y_0_xyyyz_xxxyy, g_y_0_xyyyz_xxxyz, g_y_0_xyyyz_xxxz, g_y_0_xyyyz_xxxzz, g_y_0_xyyyz_xxyy, g_y_0_xyyyz_xxyyy, g_y_0_xyyyz_xxyyz, g_y_0_xyyyz_xxyz, g_y_0_xyyyz_xxyzz, g_y_0_xyyyz_xxzz, g_y_0_xyyyz_xxzzz, g_y_0_xyyyz_xyyy, g_y_0_xyyyz_xyyyy, g_y_0_xyyyz_xyyyz, g_y_0_xyyyz_xyyz, g_y_0_xyyyz_xyyzz, g_y_0_xyyyz_xyzz, g_y_0_xyyyz_xyzzz, g_y_0_xyyyz_xzzz, g_y_0_xyyyz_xzzzz, g_y_0_xyyyz_yyyy, g_y_0_xyyyz_yyyz, g_y_0_xyyyz_yyzz, g_y_0_xyyyz_yzzz, g_y_0_xyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_xxxx[k] = -g_y_0_xyyyz_xxxx[k] * cd_x[k] + g_y_0_xyyyz_xxxxx[k];

                g_y_0_xxyyyz_xxxy[k] = -g_y_0_xyyyz_xxxy[k] * cd_x[k] + g_y_0_xyyyz_xxxxy[k];

                g_y_0_xxyyyz_xxxz[k] = -g_y_0_xyyyz_xxxz[k] * cd_x[k] + g_y_0_xyyyz_xxxxz[k];

                g_y_0_xxyyyz_xxyy[k] = -g_y_0_xyyyz_xxyy[k] * cd_x[k] + g_y_0_xyyyz_xxxyy[k];

                g_y_0_xxyyyz_xxyz[k] = -g_y_0_xyyyz_xxyz[k] * cd_x[k] + g_y_0_xyyyz_xxxyz[k];

                g_y_0_xxyyyz_xxzz[k] = -g_y_0_xyyyz_xxzz[k] * cd_x[k] + g_y_0_xyyyz_xxxzz[k];

                g_y_0_xxyyyz_xyyy[k] = -g_y_0_xyyyz_xyyy[k] * cd_x[k] + g_y_0_xyyyz_xxyyy[k];

                g_y_0_xxyyyz_xyyz[k] = -g_y_0_xyyyz_xyyz[k] * cd_x[k] + g_y_0_xyyyz_xxyyz[k];

                g_y_0_xxyyyz_xyzz[k] = -g_y_0_xyyyz_xyzz[k] * cd_x[k] + g_y_0_xyyyz_xxyzz[k];

                g_y_0_xxyyyz_xzzz[k] = -g_y_0_xyyyz_xzzz[k] * cd_x[k] + g_y_0_xyyyz_xxzzz[k];

                g_y_0_xxyyyz_yyyy[k] = -g_y_0_xyyyz_yyyy[k] * cd_x[k] + g_y_0_xyyyz_xyyyy[k];

                g_y_0_xxyyyz_yyyz[k] = -g_y_0_xyyyz_yyyz[k] * cd_x[k] + g_y_0_xyyyz_xyyyz[k];

                g_y_0_xxyyyz_yyzz[k] = -g_y_0_xyyyz_yyzz[k] * cd_x[k] + g_y_0_xyyyz_xyyzz[k];

                g_y_0_xxyyyz_yzzz[k] = -g_y_0_xyyyz_yzzz[k] * cd_x[k] + g_y_0_xyyyz_xyzzz[k];

                g_y_0_xxyyyz_zzzz[k] = -g_y_0_xyyyz_zzzz[k] * cd_x[k] + g_y_0_xyyyz_xzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 180);

            auto g_y_0_xxyyzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 181);

            auto g_y_0_xxyyzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 182);

            auto g_y_0_xxyyzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 183);

            auto g_y_0_xxyyzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 184);

            auto g_y_0_xxyyzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 185);

            auto g_y_0_xxyyzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 186);

            auto g_y_0_xxyyzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 187);

            auto g_y_0_xxyyzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 188);

            auto g_y_0_xxyyzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 189);

            auto g_y_0_xxyyzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 190);

            auto g_y_0_xxyyzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 191);

            auto g_y_0_xxyyzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 192);

            auto g_y_0_xxyyzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 193);

            auto g_y_0_xxyyzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 194);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyzz_xxxx, g_y_0_xxyyzz_xxxy, g_y_0_xxyyzz_xxxz, g_y_0_xxyyzz_xxyy, g_y_0_xxyyzz_xxyz, g_y_0_xxyyzz_xxzz, g_y_0_xxyyzz_xyyy, g_y_0_xxyyzz_xyyz, g_y_0_xxyyzz_xyzz, g_y_0_xxyyzz_xzzz, g_y_0_xxyyzz_yyyy, g_y_0_xxyyzz_yyyz, g_y_0_xxyyzz_yyzz, g_y_0_xxyyzz_yzzz, g_y_0_xxyyzz_zzzz, g_y_0_xyyzz_xxxx, g_y_0_xyyzz_xxxxx, g_y_0_xyyzz_xxxxy, g_y_0_xyyzz_xxxxz, g_y_0_xyyzz_xxxy, g_y_0_xyyzz_xxxyy, g_y_0_xyyzz_xxxyz, g_y_0_xyyzz_xxxz, g_y_0_xyyzz_xxxzz, g_y_0_xyyzz_xxyy, g_y_0_xyyzz_xxyyy, g_y_0_xyyzz_xxyyz, g_y_0_xyyzz_xxyz, g_y_0_xyyzz_xxyzz, g_y_0_xyyzz_xxzz, g_y_0_xyyzz_xxzzz, g_y_0_xyyzz_xyyy, g_y_0_xyyzz_xyyyy, g_y_0_xyyzz_xyyyz, g_y_0_xyyzz_xyyz, g_y_0_xyyzz_xyyzz, g_y_0_xyyzz_xyzz, g_y_0_xyyzz_xyzzz, g_y_0_xyyzz_xzzz, g_y_0_xyyzz_xzzzz, g_y_0_xyyzz_yyyy, g_y_0_xyyzz_yyyz, g_y_0_xyyzz_yyzz, g_y_0_xyyzz_yzzz, g_y_0_xyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_xxxx[k] = -g_y_0_xyyzz_xxxx[k] * cd_x[k] + g_y_0_xyyzz_xxxxx[k];

                g_y_0_xxyyzz_xxxy[k] = -g_y_0_xyyzz_xxxy[k] * cd_x[k] + g_y_0_xyyzz_xxxxy[k];

                g_y_0_xxyyzz_xxxz[k] = -g_y_0_xyyzz_xxxz[k] * cd_x[k] + g_y_0_xyyzz_xxxxz[k];

                g_y_0_xxyyzz_xxyy[k] = -g_y_0_xyyzz_xxyy[k] * cd_x[k] + g_y_0_xyyzz_xxxyy[k];

                g_y_0_xxyyzz_xxyz[k] = -g_y_0_xyyzz_xxyz[k] * cd_x[k] + g_y_0_xyyzz_xxxyz[k];

                g_y_0_xxyyzz_xxzz[k] = -g_y_0_xyyzz_xxzz[k] * cd_x[k] + g_y_0_xyyzz_xxxzz[k];

                g_y_0_xxyyzz_xyyy[k] = -g_y_0_xyyzz_xyyy[k] * cd_x[k] + g_y_0_xyyzz_xxyyy[k];

                g_y_0_xxyyzz_xyyz[k] = -g_y_0_xyyzz_xyyz[k] * cd_x[k] + g_y_0_xyyzz_xxyyz[k];

                g_y_0_xxyyzz_xyzz[k] = -g_y_0_xyyzz_xyzz[k] * cd_x[k] + g_y_0_xyyzz_xxyzz[k];

                g_y_0_xxyyzz_xzzz[k] = -g_y_0_xyyzz_xzzz[k] * cd_x[k] + g_y_0_xyyzz_xxzzz[k];

                g_y_0_xxyyzz_yyyy[k] = -g_y_0_xyyzz_yyyy[k] * cd_x[k] + g_y_0_xyyzz_xyyyy[k];

                g_y_0_xxyyzz_yyyz[k] = -g_y_0_xyyzz_yyyz[k] * cd_x[k] + g_y_0_xyyzz_xyyyz[k];

                g_y_0_xxyyzz_yyzz[k] = -g_y_0_xyyzz_yyzz[k] * cd_x[k] + g_y_0_xyyzz_xyyzz[k];

                g_y_0_xxyyzz_yzzz[k] = -g_y_0_xyyzz_yzzz[k] * cd_x[k] + g_y_0_xyyzz_xyzzz[k];

                g_y_0_xxyyzz_zzzz[k] = -g_y_0_xyyzz_zzzz[k] * cd_x[k] + g_y_0_xyyzz_xzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 195);

            auto g_y_0_xxyzzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 196);

            auto g_y_0_xxyzzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 197);

            auto g_y_0_xxyzzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 198);

            auto g_y_0_xxyzzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 199);

            auto g_y_0_xxyzzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 200);

            auto g_y_0_xxyzzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 201);

            auto g_y_0_xxyzzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 202);

            auto g_y_0_xxyzzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 203);

            auto g_y_0_xxyzzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 204);

            auto g_y_0_xxyzzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 205);

            auto g_y_0_xxyzzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 206);

            auto g_y_0_xxyzzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 207);

            auto g_y_0_xxyzzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 208);

            auto g_y_0_xxyzzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzzz_xxxx, g_y_0_xxyzzz_xxxy, g_y_0_xxyzzz_xxxz, g_y_0_xxyzzz_xxyy, g_y_0_xxyzzz_xxyz, g_y_0_xxyzzz_xxzz, g_y_0_xxyzzz_xyyy, g_y_0_xxyzzz_xyyz, g_y_0_xxyzzz_xyzz, g_y_0_xxyzzz_xzzz, g_y_0_xxyzzz_yyyy, g_y_0_xxyzzz_yyyz, g_y_0_xxyzzz_yyzz, g_y_0_xxyzzz_yzzz, g_y_0_xxyzzz_zzzz, g_y_0_xyzzz_xxxx, g_y_0_xyzzz_xxxxx, g_y_0_xyzzz_xxxxy, g_y_0_xyzzz_xxxxz, g_y_0_xyzzz_xxxy, g_y_0_xyzzz_xxxyy, g_y_0_xyzzz_xxxyz, g_y_0_xyzzz_xxxz, g_y_0_xyzzz_xxxzz, g_y_0_xyzzz_xxyy, g_y_0_xyzzz_xxyyy, g_y_0_xyzzz_xxyyz, g_y_0_xyzzz_xxyz, g_y_0_xyzzz_xxyzz, g_y_0_xyzzz_xxzz, g_y_0_xyzzz_xxzzz, g_y_0_xyzzz_xyyy, g_y_0_xyzzz_xyyyy, g_y_0_xyzzz_xyyyz, g_y_0_xyzzz_xyyz, g_y_0_xyzzz_xyyzz, g_y_0_xyzzz_xyzz, g_y_0_xyzzz_xyzzz, g_y_0_xyzzz_xzzz, g_y_0_xyzzz_xzzzz, g_y_0_xyzzz_yyyy, g_y_0_xyzzz_yyyz, g_y_0_xyzzz_yyzz, g_y_0_xyzzz_yzzz, g_y_0_xyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_xxxx[k] = -g_y_0_xyzzz_xxxx[k] * cd_x[k] + g_y_0_xyzzz_xxxxx[k];

                g_y_0_xxyzzz_xxxy[k] = -g_y_0_xyzzz_xxxy[k] * cd_x[k] + g_y_0_xyzzz_xxxxy[k];

                g_y_0_xxyzzz_xxxz[k] = -g_y_0_xyzzz_xxxz[k] * cd_x[k] + g_y_0_xyzzz_xxxxz[k];

                g_y_0_xxyzzz_xxyy[k] = -g_y_0_xyzzz_xxyy[k] * cd_x[k] + g_y_0_xyzzz_xxxyy[k];

                g_y_0_xxyzzz_xxyz[k] = -g_y_0_xyzzz_xxyz[k] * cd_x[k] + g_y_0_xyzzz_xxxyz[k];

                g_y_0_xxyzzz_xxzz[k] = -g_y_0_xyzzz_xxzz[k] * cd_x[k] + g_y_0_xyzzz_xxxzz[k];

                g_y_0_xxyzzz_xyyy[k] = -g_y_0_xyzzz_xyyy[k] * cd_x[k] + g_y_0_xyzzz_xxyyy[k];

                g_y_0_xxyzzz_xyyz[k] = -g_y_0_xyzzz_xyyz[k] * cd_x[k] + g_y_0_xyzzz_xxyyz[k];

                g_y_0_xxyzzz_xyzz[k] = -g_y_0_xyzzz_xyzz[k] * cd_x[k] + g_y_0_xyzzz_xxyzz[k];

                g_y_0_xxyzzz_xzzz[k] = -g_y_0_xyzzz_xzzz[k] * cd_x[k] + g_y_0_xyzzz_xxzzz[k];

                g_y_0_xxyzzz_yyyy[k] = -g_y_0_xyzzz_yyyy[k] * cd_x[k] + g_y_0_xyzzz_xyyyy[k];

                g_y_0_xxyzzz_yyyz[k] = -g_y_0_xyzzz_yyyz[k] * cd_x[k] + g_y_0_xyzzz_xyyyz[k];

                g_y_0_xxyzzz_yyzz[k] = -g_y_0_xyzzz_yyzz[k] * cd_x[k] + g_y_0_xyzzz_xyyzz[k];

                g_y_0_xxyzzz_yzzz[k] = -g_y_0_xyzzz_yzzz[k] * cd_x[k] + g_y_0_xyzzz_xyzzz[k];

                g_y_0_xxyzzz_zzzz[k] = -g_y_0_xyzzz_zzzz[k] * cd_x[k] + g_y_0_xyzzz_xzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 210);

            auto g_y_0_xxzzzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 211);

            auto g_y_0_xxzzzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 212);

            auto g_y_0_xxzzzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 213);

            auto g_y_0_xxzzzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 214);

            auto g_y_0_xxzzzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 215);

            auto g_y_0_xxzzzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 216);

            auto g_y_0_xxzzzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 217);

            auto g_y_0_xxzzzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 218);

            auto g_y_0_xxzzzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 219);

            auto g_y_0_xxzzzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 220);

            auto g_y_0_xxzzzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 221);

            auto g_y_0_xxzzzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 222);

            auto g_y_0_xxzzzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 223);

            auto g_y_0_xxzzzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 224);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzzz_xxxx, g_y_0_xxzzzz_xxxy, g_y_0_xxzzzz_xxxz, g_y_0_xxzzzz_xxyy, g_y_0_xxzzzz_xxyz, g_y_0_xxzzzz_xxzz, g_y_0_xxzzzz_xyyy, g_y_0_xxzzzz_xyyz, g_y_0_xxzzzz_xyzz, g_y_0_xxzzzz_xzzz, g_y_0_xxzzzz_yyyy, g_y_0_xxzzzz_yyyz, g_y_0_xxzzzz_yyzz, g_y_0_xxzzzz_yzzz, g_y_0_xxzzzz_zzzz, g_y_0_xzzzz_xxxx, g_y_0_xzzzz_xxxxx, g_y_0_xzzzz_xxxxy, g_y_0_xzzzz_xxxxz, g_y_0_xzzzz_xxxy, g_y_0_xzzzz_xxxyy, g_y_0_xzzzz_xxxyz, g_y_0_xzzzz_xxxz, g_y_0_xzzzz_xxxzz, g_y_0_xzzzz_xxyy, g_y_0_xzzzz_xxyyy, g_y_0_xzzzz_xxyyz, g_y_0_xzzzz_xxyz, g_y_0_xzzzz_xxyzz, g_y_0_xzzzz_xxzz, g_y_0_xzzzz_xxzzz, g_y_0_xzzzz_xyyy, g_y_0_xzzzz_xyyyy, g_y_0_xzzzz_xyyyz, g_y_0_xzzzz_xyyz, g_y_0_xzzzz_xyyzz, g_y_0_xzzzz_xyzz, g_y_0_xzzzz_xyzzz, g_y_0_xzzzz_xzzz, g_y_0_xzzzz_xzzzz, g_y_0_xzzzz_yyyy, g_y_0_xzzzz_yyyz, g_y_0_xzzzz_yyzz, g_y_0_xzzzz_yzzz, g_y_0_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_xxxx[k] = -g_y_0_xzzzz_xxxx[k] * cd_x[k] + g_y_0_xzzzz_xxxxx[k];

                g_y_0_xxzzzz_xxxy[k] = -g_y_0_xzzzz_xxxy[k] * cd_x[k] + g_y_0_xzzzz_xxxxy[k];

                g_y_0_xxzzzz_xxxz[k] = -g_y_0_xzzzz_xxxz[k] * cd_x[k] + g_y_0_xzzzz_xxxxz[k];

                g_y_0_xxzzzz_xxyy[k] = -g_y_0_xzzzz_xxyy[k] * cd_x[k] + g_y_0_xzzzz_xxxyy[k];

                g_y_0_xxzzzz_xxyz[k] = -g_y_0_xzzzz_xxyz[k] * cd_x[k] + g_y_0_xzzzz_xxxyz[k];

                g_y_0_xxzzzz_xxzz[k] = -g_y_0_xzzzz_xxzz[k] * cd_x[k] + g_y_0_xzzzz_xxxzz[k];

                g_y_0_xxzzzz_xyyy[k] = -g_y_0_xzzzz_xyyy[k] * cd_x[k] + g_y_0_xzzzz_xxyyy[k];

                g_y_0_xxzzzz_xyyz[k] = -g_y_0_xzzzz_xyyz[k] * cd_x[k] + g_y_0_xzzzz_xxyyz[k];

                g_y_0_xxzzzz_xyzz[k] = -g_y_0_xzzzz_xyzz[k] * cd_x[k] + g_y_0_xzzzz_xxyzz[k];

                g_y_0_xxzzzz_xzzz[k] = -g_y_0_xzzzz_xzzz[k] * cd_x[k] + g_y_0_xzzzz_xxzzz[k];

                g_y_0_xxzzzz_yyyy[k] = -g_y_0_xzzzz_yyyy[k] * cd_x[k] + g_y_0_xzzzz_xyyyy[k];

                g_y_0_xxzzzz_yyyz[k] = -g_y_0_xzzzz_yyyz[k] * cd_x[k] + g_y_0_xzzzz_xyyyz[k];

                g_y_0_xxzzzz_yyzz[k] = -g_y_0_xzzzz_yyzz[k] * cd_x[k] + g_y_0_xzzzz_xyyzz[k];

                g_y_0_xxzzzz_yzzz[k] = -g_y_0_xzzzz_yzzz[k] * cd_x[k] + g_y_0_xzzzz_xyzzz[k];

                g_y_0_xxzzzz_zzzz[k] = -g_y_0_xzzzz_zzzz[k] * cd_x[k] + g_y_0_xzzzz_xzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 225);

            auto g_y_0_xyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 226);

            auto g_y_0_xyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 227);

            auto g_y_0_xyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 228);

            auto g_y_0_xyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 229);

            auto g_y_0_xyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 230);

            auto g_y_0_xyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 231);

            auto g_y_0_xyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 232);

            auto g_y_0_xyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 233);

            auto g_y_0_xyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 234);

            auto g_y_0_xyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 235);

            auto g_y_0_xyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 236);

            auto g_y_0_xyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 237);

            auto g_y_0_xyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 238);

            auto g_y_0_xyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 239);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyy_xxxx, g_y_0_xyyyyy_xxxy, g_y_0_xyyyyy_xxxz, g_y_0_xyyyyy_xxyy, g_y_0_xyyyyy_xxyz, g_y_0_xyyyyy_xxzz, g_y_0_xyyyyy_xyyy, g_y_0_xyyyyy_xyyz, g_y_0_xyyyyy_xyzz, g_y_0_xyyyyy_xzzz, g_y_0_xyyyyy_yyyy, g_y_0_xyyyyy_yyyz, g_y_0_xyyyyy_yyzz, g_y_0_xyyyyy_yzzz, g_y_0_xyyyyy_zzzz, g_y_0_yyyyy_xxxx, g_y_0_yyyyy_xxxxx, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_yyyy, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_xxxx[k] = -g_y_0_yyyyy_xxxx[k] * cd_x[k] + g_y_0_yyyyy_xxxxx[k];

                g_y_0_xyyyyy_xxxy[k] = -g_y_0_yyyyy_xxxy[k] * cd_x[k] + g_y_0_yyyyy_xxxxy[k];

                g_y_0_xyyyyy_xxxz[k] = -g_y_0_yyyyy_xxxz[k] * cd_x[k] + g_y_0_yyyyy_xxxxz[k];

                g_y_0_xyyyyy_xxyy[k] = -g_y_0_yyyyy_xxyy[k] * cd_x[k] + g_y_0_yyyyy_xxxyy[k];

                g_y_0_xyyyyy_xxyz[k] = -g_y_0_yyyyy_xxyz[k] * cd_x[k] + g_y_0_yyyyy_xxxyz[k];

                g_y_0_xyyyyy_xxzz[k] = -g_y_0_yyyyy_xxzz[k] * cd_x[k] + g_y_0_yyyyy_xxxzz[k];

                g_y_0_xyyyyy_xyyy[k] = -g_y_0_yyyyy_xyyy[k] * cd_x[k] + g_y_0_yyyyy_xxyyy[k];

                g_y_0_xyyyyy_xyyz[k] = -g_y_0_yyyyy_xyyz[k] * cd_x[k] + g_y_0_yyyyy_xxyyz[k];

                g_y_0_xyyyyy_xyzz[k] = -g_y_0_yyyyy_xyzz[k] * cd_x[k] + g_y_0_yyyyy_xxyzz[k];

                g_y_0_xyyyyy_xzzz[k] = -g_y_0_yyyyy_xzzz[k] * cd_x[k] + g_y_0_yyyyy_xxzzz[k];

                g_y_0_xyyyyy_yyyy[k] = -g_y_0_yyyyy_yyyy[k] * cd_x[k] + g_y_0_yyyyy_xyyyy[k];

                g_y_0_xyyyyy_yyyz[k] = -g_y_0_yyyyy_yyyz[k] * cd_x[k] + g_y_0_yyyyy_xyyyz[k];

                g_y_0_xyyyyy_yyzz[k] = -g_y_0_yyyyy_yyzz[k] * cd_x[k] + g_y_0_yyyyy_xyyzz[k];

                g_y_0_xyyyyy_yzzz[k] = -g_y_0_yyyyy_yzzz[k] * cd_x[k] + g_y_0_yyyyy_xyzzz[k];

                g_y_0_xyyyyy_zzzz[k] = -g_y_0_yyyyy_zzzz[k] * cd_x[k] + g_y_0_yyyyy_xzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 240);

            auto g_y_0_xyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 241);

            auto g_y_0_xyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 242);

            auto g_y_0_xyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 243);

            auto g_y_0_xyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 244);

            auto g_y_0_xyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 245);

            auto g_y_0_xyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 246);

            auto g_y_0_xyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 247);

            auto g_y_0_xyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 248);

            auto g_y_0_xyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 249);

            auto g_y_0_xyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 250);

            auto g_y_0_xyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 251);

            auto g_y_0_xyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 252);

            auto g_y_0_xyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 253);

            auto g_y_0_xyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 254);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyz_xxxx, g_y_0_xyyyyz_xxxy, g_y_0_xyyyyz_xxxz, g_y_0_xyyyyz_xxyy, g_y_0_xyyyyz_xxyz, g_y_0_xyyyyz_xxzz, g_y_0_xyyyyz_xyyy, g_y_0_xyyyyz_xyyz, g_y_0_xyyyyz_xyzz, g_y_0_xyyyyz_xzzz, g_y_0_xyyyyz_yyyy, g_y_0_xyyyyz_yyyz, g_y_0_xyyyyz_yyzz, g_y_0_xyyyyz_yzzz, g_y_0_xyyyyz_zzzz, g_y_0_yyyyz_xxxx, g_y_0_yyyyz_xxxxx, g_y_0_yyyyz_xxxxy, g_y_0_yyyyz_xxxxz, g_y_0_yyyyz_xxxy, g_y_0_yyyyz_xxxyy, g_y_0_yyyyz_xxxyz, g_y_0_yyyyz_xxxz, g_y_0_yyyyz_xxxzz, g_y_0_yyyyz_xxyy, g_y_0_yyyyz_xxyyy, g_y_0_yyyyz_xxyyz, g_y_0_yyyyz_xxyz, g_y_0_yyyyz_xxyzz, g_y_0_yyyyz_xxzz, g_y_0_yyyyz_xxzzz, g_y_0_yyyyz_xyyy, g_y_0_yyyyz_xyyyy, g_y_0_yyyyz_xyyyz, g_y_0_yyyyz_xyyz, g_y_0_yyyyz_xyyzz, g_y_0_yyyyz_xyzz, g_y_0_yyyyz_xyzzz, g_y_0_yyyyz_xzzz, g_y_0_yyyyz_xzzzz, g_y_0_yyyyz_yyyy, g_y_0_yyyyz_yyyz, g_y_0_yyyyz_yyzz, g_y_0_yyyyz_yzzz, g_y_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_xxxx[k] = -g_y_0_yyyyz_xxxx[k] * cd_x[k] + g_y_0_yyyyz_xxxxx[k];

                g_y_0_xyyyyz_xxxy[k] = -g_y_0_yyyyz_xxxy[k] * cd_x[k] + g_y_0_yyyyz_xxxxy[k];

                g_y_0_xyyyyz_xxxz[k] = -g_y_0_yyyyz_xxxz[k] * cd_x[k] + g_y_0_yyyyz_xxxxz[k];

                g_y_0_xyyyyz_xxyy[k] = -g_y_0_yyyyz_xxyy[k] * cd_x[k] + g_y_0_yyyyz_xxxyy[k];

                g_y_0_xyyyyz_xxyz[k] = -g_y_0_yyyyz_xxyz[k] * cd_x[k] + g_y_0_yyyyz_xxxyz[k];

                g_y_0_xyyyyz_xxzz[k] = -g_y_0_yyyyz_xxzz[k] * cd_x[k] + g_y_0_yyyyz_xxxzz[k];

                g_y_0_xyyyyz_xyyy[k] = -g_y_0_yyyyz_xyyy[k] * cd_x[k] + g_y_0_yyyyz_xxyyy[k];

                g_y_0_xyyyyz_xyyz[k] = -g_y_0_yyyyz_xyyz[k] * cd_x[k] + g_y_0_yyyyz_xxyyz[k];

                g_y_0_xyyyyz_xyzz[k] = -g_y_0_yyyyz_xyzz[k] * cd_x[k] + g_y_0_yyyyz_xxyzz[k];

                g_y_0_xyyyyz_xzzz[k] = -g_y_0_yyyyz_xzzz[k] * cd_x[k] + g_y_0_yyyyz_xxzzz[k];

                g_y_0_xyyyyz_yyyy[k] = -g_y_0_yyyyz_yyyy[k] * cd_x[k] + g_y_0_yyyyz_xyyyy[k];

                g_y_0_xyyyyz_yyyz[k] = -g_y_0_yyyyz_yyyz[k] * cd_x[k] + g_y_0_yyyyz_xyyyz[k];

                g_y_0_xyyyyz_yyzz[k] = -g_y_0_yyyyz_yyzz[k] * cd_x[k] + g_y_0_yyyyz_xyyzz[k];

                g_y_0_xyyyyz_yzzz[k] = -g_y_0_yyyyz_yzzz[k] * cd_x[k] + g_y_0_yyyyz_xyzzz[k];

                g_y_0_xyyyyz_zzzz[k] = -g_y_0_yyyyz_zzzz[k] * cd_x[k] + g_y_0_yyyyz_xzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 255);

            auto g_y_0_xyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 256);

            auto g_y_0_xyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 257);

            auto g_y_0_xyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 258);

            auto g_y_0_xyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 259);

            auto g_y_0_xyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 260);

            auto g_y_0_xyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 261);

            auto g_y_0_xyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 262);

            auto g_y_0_xyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 263);

            auto g_y_0_xyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 264);

            auto g_y_0_xyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 265);

            auto g_y_0_xyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 266);

            auto g_y_0_xyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 267);

            auto g_y_0_xyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 268);

            auto g_y_0_xyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 269);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyzz_xxxx, g_y_0_xyyyzz_xxxy, g_y_0_xyyyzz_xxxz, g_y_0_xyyyzz_xxyy, g_y_0_xyyyzz_xxyz, g_y_0_xyyyzz_xxzz, g_y_0_xyyyzz_xyyy, g_y_0_xyyyzz_xyyz, g_y_0_xyyyzz_xyzz, g_y_0_xyyyzz_xzzz, g_y_0_xyyyzz_yyyy, g_y_0_xyyyzz_yyyz, g_y_0_xyyyzz_yyzz, g_y_0_xyyyzz_yzzz, g_y_0_xyyyzz_zzzz, g_y_0_yyyzz_xxxx, g_y_0_yyyzz_xxxxx, g_y_0_yyyzz_xxxxy, g_y_0_yyyzz_xxxxz, g_y_0_yyyzz_xxxy, g_y_0_yyyzz_xxxyy, g_y_0_yyyzz_xxxyz, g_y_0_yyyzz_xxxz, g_y_0_yyyzz_xxxzz, g_y_0_yyyzz_xxyy, g_y_0_yyyzz_xxyyy, g_y_0_yyyzz_xxyyz, g_y_0_yyyzz_xxyz, g_y_0_yyyzz_xxyzz, g_y_0_yyyzz_xxzz, g_y_0_yyyzz_xxzzz, g_y_0_yyyzz_xyyy, g_y_0_yyyzz_xyyyy, g_y_0_yyyzz_xyyyz, g_y_0_yyyzz_xyyz, g_y_0_yyyzz_xyyzz, g_y_0_yyyzz_xyzz, g_y_0_yyyzz_xyzzz, g_y_0_yyyzz_xzzz, g_y_0_yyyzz_xzzzz, g_y_0_yyyzz_yyyy, g_y_0_yyyzz_yyyz, g_y_0_yyyzz_yyzz, g_y_0_yyyzz_yzzz, g_y_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_xxxx[k] = -g_y_0_yyyzz_xxxx[k] * cd_x[k] + g_y_0_yyyzz_xxxxx[k];

                g_y_0_xyyyzz_xxxy[k] = -g_y_0_yyyzz_xxxy[k] * cd_x[k] + g_y_0_yyyzz_xxxxy[k];

                g_y_0_xyyyzz_xxxz[k] = -g_y_0_yyyzz_xxxz[k] * cd_x[k] + g_y_0_yyyzz_xxxxz[k];

                g_y_0_xyyyzz_xxyy[k] = -g_y_0_yyyzz_xxyy[k] * cd_x[k] + g_y_0_yyyzz_xxxyy[k];

                g_y_0_xyyyzz_xxyz[k] = -g_y_0_yyyzz_xxyz[k] * cd_x[k] + g_y_0_yyyzz_xxxyz[k];

                g_y_0_xyyyzz_xxzz[k] = -g_y_0_yyyzz_xxzz[k] * cd_x[k] + g_y_0_yyyzz_xxxzz[k];

                g_y_0_xyyyzz_xyyy[k] = -g_y_0_yyyzz_xyyy[k] * cd_x[k] + g_y_0_yyyzz_xxyyy[k];

                g_y_0_xyyyzz_xyyz[k] = -g_y_0_yyyzz_xyyz[k] * cd_x[k] + g_y_0_yyyzz_xxyyz[k];

                g_y_0_xyyyzz_xyzz[k] = -g_y_0_yyyzz_xyzz[k] * cd_x[k] + g_y_0_yyyzz_xxyzz[k];

                g_y_0_xyyyzz_xzzz[k] = -g_y_0_yyyzz_xzzz[k] * cd_x[k] + g_y_0_yyyzz_xxzzz[k];

                g_y_0_xyyyzz_yyyy[k] = -g_y_0_yyyzz_yyyy[k] * cd_x[k] + g_y_0_yyyzz_xyyyy[k];

                g_y_0_xyyyzz_yyyz[k] = -g_y_0_yyyzz_yyyz[k] * cd_x[k] + g_y_0_yyyzz_xyyyz[k];

                g_y_0_xyyyzz_yyzz[k] = -g_y_0_yyyzz_yyzz[k] * cd_x[k] + g_y_0_yyyzz_xyyzz[k];

                g_y_0_xyyyzz_yzzz[k] = -g_y_0_yyyzz_yzzz[k] * cd_x[k] + g_y_0_yyyzz_xyzzz[k];

                g_y_0_xyyyzz_zzzz[k] = -g_y_0_yyyzz_zzzz[k] * cd_x[k] + g_y_0_yyyzz_xzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 270);

            auto g_y_0_xyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 271);

            auto g_y_0_xyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 272);

            auto g_y_0_xyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 273);

            auto g_y_0_xyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 274);

            auto g_y_0_xyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 275);

            auto g_y_0_xyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 276);

            auto g_y_0_xyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 277);

            auto g_y_0_xyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 278);

            auto g_y_0_xyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 279);

            auto g_y_0_xyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 280);

            auto g_y_0_xyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 281);

            auto g_y_0_xyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 282);

            auto g_y_0_xyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 283);

            auto g_y_0_xyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 284);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzzz_xxxx, g_y_0_xyyzzz_xxxy, g_y_0_xyyzzz_xxxz, g_y_0_xyyzzz_xxyy, g_y_0_xyyzzz_xxyz, g_y_0_xyyzzz_xxzz, g_y_0_xyyzzz_xyyy, g_y_0_xyyzzz_xyyz, g_y_0_xyyzzz_xyzz, g_y_0_xyyzzz_xzzz, g_y_0_xyyzzz_yyyy, g_y_0_xyyzzz_yyyz, g_y_0_xyyzzz_yyzz, g_y_0_xyyzzz_yzzz, g_y_0_xyyzzz_zzzz, g_y_0_yyzzz_xxxx, g_y_0_yyzzz_xxxxx, g_y_0_yyzzz_xxxxy, g_y_0_yyzzz_xxxxz, g_y_0_yyzzz_xxxy, g_y_0_yyzzz_xxxyy, g_y_0_yyzzz_xxxyz, g_y_0_yyzzz_xxxz, g_y_0_yyzzz_xxxzz, g_y_0_yyzzz_xxyy, g_y_0_yyzzz_xxyyy, g_y_0_yyzzz_xxyyz, g_y_0_yyzzz_xxyz, g_y_0_yyzzz_xxyzz, g_y_0_yyzzz_xxzz, g_y_0_yyzzz_xxzzz, g_y_0_yyzzz_xyyy, g_y_0_yyzzz_xyyyy, g_y_0_yyzzz_xyyyz, g_y_0_yyzzz_xyyz, g_y_0_yyzzz_xyyzz, g_y_0_yyzzz_xyzz, g_y_0_yyzzz_xyzzz, g_y_0_yyzzz_xzzz, g_y_0_yyzzz_xzzzz, g_y_0_yyzzz_yyyy, g_y_0_yyzzz_yyyz, g_y_0_yyzzz_yyzz, g_y_0_yyzzz_yzzz, g_y_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_xxxx[k] = -g_y_0_yyzzz_xxxx[k] * cd_x[k] + g_y_0_yyzzz_xxxxx[k];

                g_y_0_xyyzzz_xxxy[k] = -g_y_0_yyzzz_xxxy[k] * cd_x[k] + g_y_0_yyzzz_xxxxy[k];

                g_y_0_xyyzzz_xxxz[k] = -g_y_0_yyzzz_xxxz[k] * cd_x[k] + g_y_0_yyzzz_xxxxz[k];

                g_y_0_xyyzzz_xxyy[k] = -g_y_0_yyzzz_xxyy[k] * cd_x[k] + g_y_0_yyzzz_xxxyy[k];

                g_y_0_xyyzzz_xxyz[k] = -g_y_0_yyzzz_xxyz[k] * cd_x[k] + g_y_0_yyzzz_xxxyz[k];

                g_y_0_xyyzzz_xxzz[k] = -g_y_0_yyzzz_xxzz[k] * cd_x[k] + g_y_0_yyzzz_xxxzz[k];

                g_y_0_xyyzzz_xyyy[k] = -g_y_0_yyzzz_xyyy[k] * cd_x[k] + g_y_0_yyzzz_xxyyy[k];

                g_y_0_xyyzzz_xyyz[k] = -g_y_0_yyzzz_xyyz[k] * cd_x[k] + g_y_0_yyzzz_xxyyz[k];

                g_y_0_xyyzzz_xyzz[k] = -g_y_0_yyzzz_xyzz[k] * cd_x[k] + g_y_0_yyzzz_xxyzz[k];

                g_y_0_xyyzzz_xzzz[k] = -g_y_0_yyzzz_xzzz[k] * cd_x[k] + g_y_0_yyzzz_xxzzz[k];

                g_y_0_xyyzzz_yyyy[k] = -g_y_0_yyzzz_yyyy[k] * cd_x[k] + g_y_0_yyzzz_xyyyy[k];

                g_y_0_xyyzzz_yyyz[k] = -g_y_0_yyzzz_yyyz[k] * cd_x[k] + g_y_0_yyzzz_xyyyz[k];

                g_y_0_xyyzzz_yyzz[k] = -g_y_0_yyzzz_yyzz[k] * cd_x[k] + g_y_0_yyzzz_xyyzz[k];

                g_y_0_xyyzzz_yzzz[k] = -g_y_0_yyzzz_yzzz[k] * cd_x[k] + g_y_0_yyzzz_xyzzz[k];

                g_y_0_xyyzzz_zzzz[k] = -g_y_0_yyzzz_zzzz[k] * cd_x[k] + g_y_0_yyzzz_xzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 285);

            auto g_y_0_xyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 286);

            auto g_y_0_xyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 287);

            auto g_y_0_xyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 288);

            auto g_y_0_xyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 289);

            auto g_y_0_xyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 290);

            auto g_y_0_xyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 291);

            auto g_y_0_xyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 292);

            auto g_y_0_xyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 293);

            auto g_y_0_xyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 294);

            auto g_y_0_xyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 295);

            auto g_y_0_xyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 296);

            auto g_y_0_xyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 297);

            auto g_y_0_xyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 298);

            auto g_y_0_xyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 299);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzzz_xxxx, g_y_0_xyzzzz_xxxy, g_y_0_xyzzzz_xxxz, g_y_0_xyzzzz_xxyy, g_y_0_xyzzzz_xxyz, g_y_0_xyzzzz_xxzz, g_y_0_xyzzzz_xyyy, g_y_0_xyzzzz_xyyz, g_y_0_xyzzzz_xyzz, g_y_0_xyzzzz_xzzz, g_y_0_xyzzzz_yyyy, g_y_0_xyzzzz_yyyz, g_y_0_xyzzzz_yyzz, g_y_0_xyzzzz_yzzz, g_y_0_xyzzzz_zzzz, g_y_0_yzzzz_xxxx, g_y_0_yzzzz_xxxxx, g_y_0_yzzzz_xxxxy, g_y_0_yzzzz_xxxxz, g_y_0_yzzzz_xxxy, g_y_0_yzzzz_xxxyy, g_y_0_yzzzz_xxxyz, g_y_0_yzzzz_xxxz, g_y_0_yzzzz_xxxzz, g_y_0_yzzzz_xxyy, g_y_0_yzzzz_xxyyy, g_y_0_yzzzz_xxyyz, g_y_0_yzzzz_xxyz, g_y_0_yzzzz_xxyzz, g_y_0_yzzzz_xxzz, g_y_0_yzzzz_xxzzz, g_y_0_yzzzz_xyyy, g_y_0_yzzzz_xyyyy, g_y_0_yzzzz_xyyyz, g_y_0_yzzzz_xyyz, g_y_0_yzzzz_xyyzz, g_y_0_yzzzz_xyzz, g_y_0_yzzzz_xyzzz, g_y_0_yzzzz_xzzz, g_y_0_yzzzz_xzzzz, g_y_0_yzzzz_yyyy, g_y_0_yzzzz_yyyz, g_y_0_yzzzz_yyzz, g_y_0_yzzzz_yzzz, g_y_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_xxxx[k] = -g_y_0_yzzzz_xxxx[k] * cd_x[k] + g_y_0_yzzzz_xxxxx[k];

                g_y_0_xyzzzz_xxxy[k] = -g_y_0_yzzzz_xxxy[k] * cd_x[k] + g_y_0_yzzzz_xxxxy[k];

                g_y_0_xyzzzz_xxxz[k] = -g_y_0_yzzzz_xxxz[k] * cd_x[k] + g_y_0_yzzzz_xxxxz[k];

                g_y_0_xyzzzz_xxyy[k] = -g_y_0_yzzzz_xxyy[k] * cd_x[k] + g_y_0_yzzzz_xxxyy[k];

                g_y_0_xyzzzz_xxyz[k] = -g_y_0_yzzzz_xxyz[k] * cd_x[k] + g_y_0_yzzzz_xxxyz[k];

                g_y_0_xyzzzz_xxzz[k] = -g_y_0_yzzzz_xxzz[k] * cd_x[k] + g_y_0_yzzzz_xxxzz[k];

                g_y_0_xyzzzz_xyyy[k] = -g_y_0_yzzzz_xyyy[k] * cd_x[k] + g_y_0_yzzzz_xxyyy[k];

                g_y_0_xyzzzz_xyyz[k] = -g_y_0_yzzzz_xyyz[k] * cd_x[k] + g_y_0_yzzzz_xxyyz[k];

                g_y_0_xyzzzz_xyzz[k] = -g_y_0_yzzzz_xyzz[k] * cd_x[k] + g_y_0_yzzzz_xxyzz[k];

                g_y_0_xyzzzz_xzzz[k] = -g_y_0_yzzzz_xzzz[k] * cd_x[k] + g_y_0_yzzzz_xxzzz[k];

                g_y_0_xyzzzz_yyyy[k] = -g_y_0_yzzzz_yyyy[k] * cd_x[k] + g_y_0_yzzzz_xyyyy[k];

                g_y_0_xyzzzz_yyyz[k] = -g_y_0_yzzzz_yyyz[k] * cd_x[k] + g_y_0_yzzzz_xyyyz[k];

                g_y_0_xyzzzz_yyzz[k] = -g_y_0_yzzzz_yyzz[k] * cd_x[k] + g_y_0_yzzzz_xyyzz[k];

                g_y_0_xyzzzz_yzzz[k] = -g_y_0_yzzzz_yzzz[k] * cd_x[k] + g_y_0_yzzzz_xyzzz[k];

                g_y_0_xyzzzz_zzzz[k] = -g_y_0_yzzzz_zzzz[k] * cd_x[k] + g_y_0_yzzzz_xzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 300);

            auto g_y_0_xzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 301);

            auto g_y_0_xzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 302);

            auto g_y_0_xzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 303);

            auto g_y_0_xzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 304);

            auto g_y_0_xzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 305);

            auto g_y_0_xzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 306);

            auto g_y_0_xzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 307);

            auto g_y_0_xzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 308);

            auto g_y_0_xzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 309);

            auto g_y_0_xzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 310);

            auto g_y_0_xzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 311);

            auto g_y_0_xzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 312);

            auto g_y_0_xzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 313);

            auto g_y_0_xzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzzz_xxxx, g_y_0_xzzzzz_xxxy, g_y_0_xzzzzz_xxxz, g_y_0_xzzzzz_xxyy, g_y_0_xzzzzz_xxyz, g_y_0_xzzzzz_xxzz, g_y_0_xzzzzz_xyyy, g_y_0_xzzzzz_xyyz, g_y_0_xzzzzz_xyzz, g_y_0_xzzzzz_xzzz, g_y_0_xzzzzz_yyyy, g_y_0_xzzzzz_yyyz, g_y_0_xzzzzz_yyzz, g_y_0_xzzzzz_yzzz, g_y_0_xzzzzz_zzzz, g_y_0_zzzzz_xxxx, g_y_0_zzzzz_xxxxx, g_y_0_zzzzz_xxxxy, g_y_0_zzzzz_xxxxz, g_y_0_zzzzz_xxxy, g_y_0_zzzzz_xxxyy, g_y_0_zzzzz_xxxyz, g_y_0_zzzzz_xxxz, g_y_0_zzzzz_xxxzz, g_y_0_zzzzz_xxyy, g_y_0_zzzzz_xxyyy, g_y_0_zzzzz_xxyyz, g_y_0_zzzzz_xxyz, g_y_0_zzzzz_xxyzz, g_y_0_zzzzz_xxzz, g_y_0_zzzzz_xxzzz, g_y_0_zzzzz_xyyy, g_y_0_zzzzz_xyyyy, g_y_0_zzzzz_xyyyz, g_y_0_zzzzz_xyyz, g_y_0_zzzzz_xyyzz, g_y_0_zzzzz_xyzz, g_y_0_zzzzz_xyzzz, g_y_0_zzzzz_xzzz, g_y_0_zzzzz_xzzzz, g_y_0_zzzzz_yyyy, g_y_0_zzzzz_yyyz, g_y_0_zzzzz_yyzz, g_y_0_zzzzz_yzzz, g_y_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_xxxx[k] = -g_y_0_zzzzz_xxxx[k] * cd_x[k] + g_y_0_zzzzz_xxxxx[k];

                g_y_0_xzzzzz_xxxy[k] = -g_y_0_zzzzz_xxxy[k] * cd_x[k] + g_y_0_zzzzz_xxxxy[k];

                g_y_0_xzzzzz_xxxz[k] = -g_y_0_zzzzz_xxxz[k] * cd_x[k] + g_y_0_zzzzz_xxxxz[k];

                g_y_0_xzzzzz_xxyy[k] = -g_y_0_zzzzz_xxyy[k] * cd_x[k] + g_y_0_zzzzz_xxxyy[k];

                g_y_0_xzzzzz_xxyz[k] = -g_y_0_zzzzz_xxyz[k] * cd_x[k] + g_y_0_zzzzz_xxxyz[k];

                g_y_0_xzzzzz_xxzz[k] = -g_y_0_zzzzz_xxzz[k] * cd_x[k] + g_y_0_zzzzz_xxxzz[k];

                g_y_0_xzzzzz_xyyy[k] = -g_y_0_zzzzz_xyyy[k] * cd_x[k] + g_y_0_zzzzz_xxyyy[k];

                g_y_0_xzzzzz_xyyz[k] = -g_y_0_zzzzz_xyyz[k] * cd_x[k] + g_y_0_zzzzz_xxyyz[k];

                g_y_0_xzzzzz_xyzz[k] = -g_y_0_zzzzz_xyzz[k] * cd_x[k] + g_y_0_zzzzz_xxyzz[k];

                g_y_0_xzzzzz_xzzz[k] = -g_y_0_zzzzz_xzzz[k] * cd_x[k] + g_y_0_zzzzz_xxzzz[k];

                g_y_0_xzzzzz_yyyy[k] = -g_y_0_zzzzz_yyyy[k] * cd_x[k] + g_y_0_zzzzz_xyyyy[k];

                g_y_0_xzzzzz_yyyz[k] = -g_y_0_zzzzz_yyyz[k] * cd_x[k] + g_y_0_zzzzz_xyyyz[k];

                g_y_0_xzzzzz_yyzz[k] = -g_y_0_zzzzz_yyzz[k] * cd_x[k] + g_y_0_zzzzz_xyyzz[k];

                g_y_0_xzzzzz_yzzz[k] = -g_y_0_zzzzz_yzzz[k] * cd_x[k] + g_y_0_zzzzz_xyzzz[k];

                g_y_0_xzzzzz_zzzz[k] = -g_y_0_zzzzz_zzzz[k] * cd_x[k] + g_y_0_zzzzz_xzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 315);

            auto g_y_0_yyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 316);

            auto g_y_0_yyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 317);

            auto g_y_0_yyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 318);

            auto g_y_0_yyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 319);

            auto g_y_0_yyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 320);

            auto g_y_0_yyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 321);

            auto g_y_0_yyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 322);

            auto g_y_0_yyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 323);

            auto g_y_0_yyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 324);

            auto g_y_0_yyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 325);

            auto g_y_0_yyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 326);

            auto g_y_0_yyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 327);

            auto g_y_0_yyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 328);

            auto g_y_0_yyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 329);

            #pragma omp simd aligned(cd_y, g_y_0_yyyyy_xxxx, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_yyyy, g_y_0_yyyyy_yyyyy, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_zzzz, g_y_0_yyyyyy_xxxx, g_y_0_yyyyyy_xxxy, g_y_0_yyyyyy_xxxz, g_y_0_yyyyyy_xxyy, g_y_0_yyyyyy_xxyz, g_y_0_yyyyyy_xxzz, g_y_0_yyyyyy_xyyy, g_y_0_yyyyyy_xyyz, g_y_0_yyyyyy_xyzz, g_y_0_yyyyyy_xzzz, g_y_0_yyyyyy_yyyy, g_y_0_yyyyyy_yyyz, g_y_0_yyyyyy_yyzz, g_y_0_yyyyyy_yzzz, g_y_0_yyyyyy_zzzz, g_yyyyy_xxxx, g_yyyyy_xxxy, g_yyyyy_xxxz, g_yyyyy_xxyy, g_yyyyy_xxyz, g_yyyyy_xxzz, g_yyyyy_xyyy, g_yyyyy_xyyz, g_yyyyy_xyzz, g_yyyyy_xzzz, g_yyyyy_yyyy, g_yyyyy_yyyz, g_yyyyy_yyzz, g_yyyyy_yzzz, g_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_xxxx[k] = -g_yyyyy_xxxx[k] - g_y_0_yyyyy_xxxx[k] * cd_y[k] + g_y_0_yyyyy_xxxxy[k];

                g_y_0_yyyyyy_xxxy[k] = -g_yyyyy_xxxy[k] - g_y_0_yyyyy_xxxy[k] * cd_y[k] + g_y_0_yyyyy_xxxyy[k];

                g_y_0_yyyyyy_xxxz[k] = -g_yyyyy_xxxz[k] - g_y_0_yyyyy_xxxz[k] * cd_y[k] + g_y_0_yyyyy_xxxyz[k];

                g_y_0_yyyyyy_xxyy[k] = -g_yyyyy_xxyy[k] - g_y_0_yyyyy_xxyy[k] * cd_y[k] + g_y_0_yyyyy_xxyyy[k];

                g_y_0_yyyyyy_xxyz[k] = -g_yyyyy_xxyz[k] - g_y_0_yyyyy_xxyz[k] * cd_y[k] + g_y_0_yyyyy_xxyyz[k];

                g_y_0_yyyyyy_xxzz[k] = -g_yyyyy_xxzz[k] - g_y_0_yyyyy_xxzz[k] * cd_y[k] + g_y_0_yyyyy_xxyzz[k];

                g_y_0_yyyyyy_xyyy[k] = -g_yyyyy_xyyy[k] - g_y_0_yyyyy_xyyy[k] * cd_y[k] + g_y_0_yyyyy_xyyyy[k];

                g_y_0_yyyyyy_xyyz[k] = -g_yyyyy_xyyz[k] - g_y_0_yyyyy_xyyz[k] * cd_y[k] + g_y_0_yyyyy_xyyyz[k];

                g_y_0_yyyyyy_xyzz[k] = -g_yyyyy_xyzz[k] - g_y_0_yyyyy_xyzz[k] * cd_y[k] + g_y_0_yyyyy_xyyzz[k];

                g_y_0_yyyyyy_xzzz[k] = -g_yyyyy_xzzz[k] - g_y_0_yyyyy_xzzz[k] * cd_y[k] + g_y_0_yyyyy_xyzzz[k];

                g_y_0_yyyyyy_yyyy[k] = -g_yyyyy_yyyy[k] - g_y_0_yyyyy_yyyy[k] * cd_y[k] + g_y_0_yyyyy_yyyyy[k];

                g_y_0_yyyyyy_yyyz[k] = -g_yyyyy_yyyz[k] - g_y_0_yyyyy_yyyz[k] * cd_y[k] + g_y_0_yyyyy_yyyyz[k];

                g_y_0_yyyyyy_yyzz[k] = -g_yyyyy_yyzz[k] - g_y_0_yyyyy_yyzz[k] * cd_y[k] + g_y_0_yyyyy_yyyzz[k];

                g_y_0_yyyyyy_yzzz[k] = -g_yyyyy_yzzz[k] - g_y_0_yyyyy_yzzz[k] * cd_y[k] + g_y_0_yyyyy_yyzzz[k];

                g_y_0_yyyyyy_zzzz[k] = -g_yyyyy_zzzz[k] - g_y_0_yyyyy_zzzz[k] * cd_y[k] + g_y_0_yyyyy_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 330);

            auto g_y_0_yyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 331);

            auto g_y_0_yyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 332);

            auto g_y_0_yyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 333);

            auto g_y_0_yyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 334);

            auto g_y_0_yyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 335);

            auto g_y_0_yyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 336);

            auto g_y_0_yyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 337);

            auto g_y_0_yyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 338);

            auto g_y_0_yyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 339);

            auto g_y_0_yyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 340);

            auto g_y_0_yyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 341);

            auto g_y_0_yyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 342);

            auto g_y_0_yyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 343);

            auto g_y_0_yyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 344);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyy_xxxx, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_yyyy, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_zzzz, g_y_0_yyyyy_zzzzz, g_y_0_yyyyyz_xxxx, g_y_0_yyyyyz_xxxy, g_y_0_yyyyyz_xxxz, g_y_0_yyyyyz_xxyy, g_y_0_yyyyyz_xxyz, g_y_0_yyyyyz_xxzz, g_y_0_yyyyyz_xyyy, g_y_0_yyyyyz_xyyz, g_y_0_yyyyyz_xyzz, g_y_0_yyyyyz_xzzz, g_y_0_yyyyyz_yyyy, g_y_0_yyyyyz_yyyz, g_y_0_yyyyyz_yyzz, g_y_0_yyyyyz_yzzz, g_y_0_yyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_xxxx[k] = -g_y_0_yyyyy_xxxx[k] * cd_z[k] + g_y_0_yyyyy_xxxxz[k];

                g_y_0_yyyyyz_xxxy[k] = -g_y_0_yyyyy_xxxy[k] * cd_z[k] + g_y_0_yyyyy_xxxyz[k];

                g_y_0_yyyyyz_xxxz[k] = -g_y_0_yyyyy_xxxz[k] * cd_z[k] + g_y_0_yyyyy_xxxzz[k];

                g_y_0_yyyyyz_xxyy[k] = -g_y_0_yyyyy_xxyy[k] * cd_z[k] + g_y_0_yyyyy_xxyyz[k];

                g_y_0_yyyyyz_xxyz[k] = -g_y_0_yyyyy_xxyz[k] * cd_z[k] + g_y_0_yyyyy_xxyzz[k];

                g_y_0_yyyyyz_xxzz[k] = -g_y_0_yyyyy_xxzz[k] * cd_z[k] + g_y_0_yyyyy_xxzzz[k];

                g_y_0_yyyyyz_xyyy[k] = -g_y_0_yyyyy_xyyy[k] * cd_z[k] + g_y_0_yyyyy_xyyyz[k];

                g_y_0_yyyyyz_xyyz[k] = -g_y_0_yyyyy_xyyz[k] * cd_z[k] + g_y_0_yyyyy_xyyzz[k];

                g_y_0_yyyyyz_xyzz[k] = -g_y_0_yyyyy_xyzz[k] * cd_z[k] + g_y_0_yyyyy_xyzzz[k];

                g_y_0_yyyyyz_xzzz[k] = -g_y_0_yyyyy_xzzz[k] * cd_z[k] + g_y_0_yyyyy_xzzzz[k];

                g_y_0_yyyyyz_yyyy[k] = -g_y_0_yyyyy_yyyy[k] * cd_z[k] + g_y_0_yyyyy_yyyyz[k];

                g_y_0_yyyyyz_yyyz[k] = -g_y_0_yyyyy_yyyz[k] * cd_z[k] + g_y_0_yyyyy_yyyzz[k];

                g_y_0_yyyyyz_yyzz[k] = -g_y_0_yyyyy_yyzz[k] * cd_z[k] + g_y_0_yyyyy_yyzzz[k];

                g_y_0_yyyyyz_yzzz[k] = -g_y_0_yyyyy_yzzz[k] * cd_z[k] + g_y_0_yyyyy_yzzzz[k];

                g_y_0_yyyyyz_zzzz[k] = -g_y_0_yyyyy_zzzz[k] * cd_z[k] + g_y_0_yyyyy_zzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 345);

            auto g_y_0_yyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 346);

            auto g_y_0_yyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 347);

            auto g_y_0_yyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 348);

            auto g_y_0_yyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 349);

            auto g_y_0_yyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 350);

            auto g_y_0_yyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 351);

            auto g_y_0_yyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 352);

            auto g_y_0_yyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 353);

            auto g_y_0_yyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 354);

            auto g_y_0_yyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 355);

            auto g_y_0_yyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 356);

            auto g_y_0_yyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 357);

            auto g_y_0_yyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 358);

            auto g_y_0_yyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 359);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyz_xxxx, g_y_0_yyyyz_xxxxz, g_y_0_yyyyz_xxxy, g_y_0_yyyyz_xxxyz, g_y_0_yyyyz_xxxz, g_y_0_yyyyz_xxxzz, g_y_0_yyyyz_xxyy, g_y_0_yyyyz_xxyyz, g_y_0_yyyyz_xxyz, g_y_0_yyyyz_xxyzz, g_y_0_yyyyz_xxzz, g_y_0_yyyyz_xxzzz, g_y_0_yyyyz_xyyy, g_y_0_yyyyz_xyyyz, g_y_0_yyyyz_xyyz, g_y_0_yyyyz_xyyzz, g_y_0_yyyyz_xyzz, g_y_0_yyyyz_xyzzz, g_y_0_yyyyz_xzzz, g_y_0_yyyyz_xzzzz, g_y_0_yyyyz_yyyy, g_y_0_yyyyz_yyyyz, g_y_0_yyyyz_yyyz, g_y_0_yyyyz_yyyzz, g_y_0_yyyyz_yyzz, g_y_0_yyyyz_yyzzz, g_y_0_yyyyz_yzzz, g_y_0_yyyyz_yzzzz, g_y_0_yyyyz_zzzz, g_y_0_yyyyz_zzzzz, g_y_0_yyyyzz_xxxx, g_y_0_yyyyzz_xxxy, g_y_0_yyyyzz_xxxz, g_y_0_yyyyzz_xxyy, g_y_0_yyyyzz_xxyz, g_y_0_yyyyzz_xxzz, g_y_0_yyyyzz_xyyy, g_y_0_yyyyzz_xyyz, g_y_0_yyyyzz_xyzz, g_y_0_yyyyzz_xzzz, g_y_0_yyyyzz_yyyy, g_y_0_yyyyzz_yyyz, g_y_0_yyyyzz_yyzz, g_y_0_yyyyzz_yzzz, g_y_0_yyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_xxxx[k] = -g_y_0_yyyyz_xxxx[k] * cd_z[k] + g_y_0_yyyyz_xxxxz[k];

                g_y_0_yyyyzz_xxxy[k] = -g_y_0_yyyyz_xxxy[k] * cd_z[k] + g_y_0_yyyyz_xxxyz[k];

                g_y_0_yyyyzz_xxxz[k] = -g_y_0_yyyyz_xxxz[k] * cd_z[k] + g_y_0_yyyyz_xxxzz[k];

                g_y_0_yyyyzz_xxyy[k] = -g_y_0_yyyyz_xxyy[k] * cd_z[k] + g_y_0_yyyyz_xxyyz[k];

                g_y_0_yyyyzz_xxyz[k] = -g_y_0_yyyyz_xxyz[k] * cd_z[k] + g_y_0_yyyyz_xxyzz[k];

                g_y_0_yyyyzz_xxzz[k] = -g_y_0_yyyyz_xxzz[k] * cd_z[k] + g_y_0_yyyyz_xxzzz[k];

                g_y_0_yyyyzz_xyyy[k] = -g_y_0_yyyyz_xyyy[k] * cd_z[k] + g_y_0_yyyyz_xyyyz[k];

                g_y_0_yyyyzz_xyyz[k] = -g_y_0_yyyyz_xyyz[k] * cd_z[k] + g_y_0_yyyyz_xyyzz[k];

                g_y_0_yyyyzz_xyzz[k] = -g_y_0_yyyyz_xyzz[k] * cd_z[k] + g_y_0_yyyyz_xyzzz[k];

                g_y_0_yyyyzz_xzzz[k] = -g_y_0_yyyyz_xzzz[k] * cd_z[k] + g_y_0_yyyyz_xzzzz[k];

                g_y_0_yyyyzz_yyyy[k] = -g_y_0_yyyyz_yyyy[k] * cd_z[k] + g_y_0_yyyyz_yyyyz[k];

                g_y_0_yyyyzz_yyyz[k] = -g_y_0_yyyyz_yyyz[k] * cd_z[k] + g_y_0_yyyyz_yyyzz[k];

                g_y_0_yyyyzz_yyzz[k] = -g_y_0_yyyyz_yyzz[k] * cd_z[k] + g_y_0_yyyyz_yyzzz[k];

                g_y_0_yyyyzz_yzzz[k] = -g_y_0_yyyyz_yzzz[k] * cd_z[k] + g_y_0_yyyyz_yzzzz[k];

                g_y_0_yyyyzz_zzzz[k] = -g_y_0_yyyyz_zzzz[k] * cd_z[k] + g_y_0_yyyyz_zzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 360);

            auto g_y_0_yyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 361);

            auto g_y_0_yyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 362);

            auto g_y_0_yyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 363);

            auto g_y_0_yyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 364);

            auto g_y_0_yyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 365);

            auto g_y_0_yyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 366);

            auto g_y_0_yyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 367);

            auto g_y_0_yyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 368);

            auto g_y_0_yyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 369);

            auto g_y_0_yyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 370);

            auto g_y_0_yyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 371);

            auto g_y_0_yyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 372);

            auto g_y_0_yyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 373);

            auto g_y_0_yyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 374);

            #pragma omp simd aligned(cd_z, g_y_0_yyyzz_xxxx, g_y_0_yyyzz_xxxxz, g_y_0_yyyzz_xxxy, g_y_0_yyyzz_xxxyz, g_y_0_yyyzz_xxxz, g_y_0_yyyzz_xxxzz, g_y_0_yyyzz_xxyy, g_y_0_yyyzz_xxyyz, g_y_0_yyyzz_xxyz, g_y_0_yyyzz_xxyzz, g_y_0_yyyzz_xxzz, g_y_0_yyyzz_xxzzz, g_y_0_yyyzz_xyyy, g_y_0_yyyzz_xyyyz, g_y_0_yyyzz_xyyz, g_y_0_yyyzz_xyyzz, g_y_0_yyyzz_xyzz, g_y_0_yyyzz_xyzzz, g_y_0_yyyzz_xzzz, g_y_0_yyyzz_xzzzz, g_y_0_yyyzz_yyyy, g_y_0_yyyzz_yyyyz, g_y_0_yyyzz_yyyz, g_y_0_yyyzz_yyyzz, g_y_0_yyyzz_yyzz, g_y_0_yyyzz_yyzzz, g_y_0_yyyzz_yzzz, g_y_0_yyyzz_yzzzz, g_y_0_yyyzz_zzzz, g_y_0_yyyzz_zzzzz, g_y_0_yyyzzz_xxxx, g_y_0_yyyzzz_xxxy, g_y_0_yyyzzz_xxxz, g_y_0_yyyzzz_xxyy, g_y_0_yyyzzz_xxyz, g_y_0_yyyzzz_xxzz, g_y_0_yyyzzz_xyyy, g_y_0_yyyzzz_xyyz, g_y_0_yyyzzz_xyzz, g_y_0_yyyzzz_xzzz, g_y_0_yyyzzz_yyyy, g_y_0_yyyzzz_yyyz, g_y_0_yyyzzz_yyzz, g_y_0_yyyzzz_yzzz, g_y_0_yyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_xxxx[k] = -g_y_0_yyyzz_xxxx[k] * cd_z[k] + g_y_0_yyyzz_xxxxz[k];

                g_y_0_yyyzzz_xxxy[k] = -g_y_0_yyyzz_xxxy[k] * cd_z[k] + g_y_0_yyyzz_xxxyz[k];

                g_y_0_yyyzzz_xxxz[k] = -g_y_0_yyyzz_xxxz[k] * cd_z[k] + g_y_0_yyyzz_xxxzz[k];

                g_y_0_yyyzzz_xxyy[k] = -g_y_0_yyyzz_xxyy[k] * cd_z[k] + g_y_0_yyyzz_xxyyz[k];

                g_y_0_yyyzzz_xxyz[k] = -g_y_0_yyyzz_xxyz[k] * cd_z[k] + g_y_0_yyyzz_xxyzz[k];

                g_y_0_yyyzzz_xxzz[k] = -g_y_0_yyyzz_xxzz[k] * cd_z[k] + g_y_0_yyyzz_xxzzz[k];

                g_y_0_yyyzzz_xyyy[k] = -g_y_0_yyyzz_xyyy[k] * cd_z[k] + g_y_0_yyyzz_xyyyz[k];

                g_y_0_yyyzzz_xyyz[k] = -g_y_0_yyyzz_xyyz[k] * cd_z[k] + g_y_0_yyyzz_xyyzz[k];

                g_y_0_yyyzzz_xyzz[k] = -g_y_0_yyyzz_xyzz[k] * cd_z[k] + g_y_0_yyyzz_xyzzz[k];

                g_y_0_yyyzzz_xzzz[k] = -g_y_0_yyyzz_xzzz[k] * cd_z[k] + g_y_0_yyyzz_xzzzz[k];

                g_y_0_yyyzzz_yyyy[k] = -g_y_0_yyyzz_yyyy[k] * cd_z[k] + g_y_0_yyyzz_yyyyz[k];

                g_y_0_yyyzzz_yyyz[k] = -g_y_0_yyyzz_yyyz[k] * cd_z[k] + g_y_0_yyyzz_yyyzz[k];

                g_y_0_yyyzzz_yyzz[k] = -g_y_0_yyyzz_yyzz[k] * cd_z[k] + g_y_0_yyyzz_yyzzz[k];

                g_y_0_yyyzzz_yzzz[k] = -g_y_0_yyyzz_yzzz[k] * cd_z[k] + g_y_0_yyyzz_yzzzz[k];

                g_y_0_yyyzzz_zzzz[k] = -g_y_0_yyyzz_zzzz[k] * cd_z[k] + g_y_0_yyyzz_zzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 375);

            auto g_y_0_yyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 376);

            auto g_y_0_yyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 377);

            auto g_y_0_yyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 378);

            auto g_y_0_yyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 379);

            auto g_y_0_yyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 380);

            auto g_y_0_yyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 381);

            auto g_y_0_yyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 382);

            auto g_y_0_yyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 383);

            auto g_y_0_yyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 384);

            auto g_y_0_yyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 385);

            auto g_y_0_yyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 386);

            auto g_y_0_yyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 387);

            auto g_y_0_yyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 388);

            auto g_y_0_yyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 389);

            #pragma omp simd aligned(cd_z, g_y_0_yyzzz_xxxx, g_y_0_yyzzz_xxxxz, g_y_0_yyzzz_xxxy, g_y_0_yyzzz_xxxyz, g_y_0_yyzzz_xxxz, g_y_0_yyzzz_xxxzz, g_y_0_yyzzz_xxyy, g_y_0_yyzzz_xxyyz, g_y_0_yyzzz_xxyz, g_y_0_yyzzz_xxyzz, g_y_0_yyzzz_xxzz, g_y_0_yyzzz_xxzzz, g_y_0_yyzzz_xyyy, g_y_0_yyzzz_xyyyz, g_y_0_yyzzz_xyyz, g_y_0_yyzzz_xyyzz, g_y_0_yyzzz_xyzz, g_y_0_yyzzz_xyzzz, g_y_0_yyzzz_xzzz, g_y_0_yyzzz_xzzzz, g_y_0_yyzzz_yyyy, g_y_0_yyzzz_yyyyz, g_y_0_yyzzz_yyyz, g_y_0_yyzzz_yyyzz, g_y_0_yyzzz_yyzz, g_y_0_yyzzz_yyzzz, g_y_0_yyzzz_yzzz, g_y_0_yyzzz_yzzzz, g_y_0_yyzzz_zzzz, g_y_0_yyzzz_zzzzz, g_y_0_yyzzzz_xxxx, g_y_0_yyzzzz_xxxy, g_y_0_yyzzzz_xxxz, g_y_0_yyzzzz_xxyy, g_y_0_yyzzzz_xxyz, g_y_0_yyzzzz_xxzz, g_y_0_yyzzzz_xyyy, g_y_0_yyzzzz_xyyz, g_y_0_yyzzzz_xyzz, g_y_0_yyzzzz_xzzz, g_y_0_yyzzzz_yyyy, g_y_0_yyzzzz_yyyz, g_y_0_yyzzzz_yyzz, g_y_0_yyzzzz_yzzz, g_y_0_yyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_xxxx[k] = -g_y_0_yyzzz_xxxx[k] * cd_z[k] + g_y_0_yyzzz_xxxxz[k];

                g_y_0_yyzzzz_xxxy[k] = -g_y_0_yyzzz_xxxy[k] * cd_z[k] + g_y_0_yyzzz_xxxyz[k];

                g_y_0_yyzzzz_xxxz[k] = -g_y_0_yyzzz_xxxz[k] * cd_z[k] + g_y_0_yyzzz_xxxzz[k];

                g_y_0_yyzzzz_xxyy[k] = -g_y_0_yyzzz_xxyy[k] * cd_z[k] + g_y_0_yyzzz_xxyyz[k];

                g_y_0_yyzzzz_xxyz[k] = -g_y_0_yyzzz_xxyz[k] * cd_z[k] + g_y_0_yyzzz_xxyzz[k];

                g_y_0_yyzzzz_xxzz[k] = -g_y_0_yyzzz_xxzz[k] * cd_z[k] + g_y_0_yyzzz_xxzzz[k];

                g_y_0_yyzzzz_xyyy[k] = -g_y_0_yyzzz_xyyy[k] * cd_z[k] + g_y_0_yyzzz_xyyyz[k];

                g_y_0_yyzzzz_xyyz[k] = -g_y_0_yyzzz_xyyz[k] * cd_z[k] + g_y_0_yyzzz_xyyzz[k];

                g_y_0_yyzzzz_xyzz[k] = -g_y_0_yyzzz_xyzz[k] * cd_z[k] + g_y_0_yyzzz_xyzzz[k];

                g_y_0_yyzzzz_xzzz[k] = -g_y_0_yyzzz_xzzz[k] * cd_z[k] + g_y_0_yyzzz_xzzzz[k];

                g_y_0_yyzzzz_yyyy[k] = -g_y_0_yyzzz_yyyy[k] * cd_z[k] + g_y_0_yyzzz_yyyyz[k];

                g_y_0_yyzzzz_yyyz[k] = -g_y_0_yyzzz_yyyz[k] * cd_z[k] + g_y_0_yyzzz_yyyzz[k];

                g_y_0_yyzzzz_yyzz[k] = -g_y_0_yyzzz_yyzz[k] * cd_z[k] + g_y_0_yyzzz_yyzzz[k];

                g_y_0_yyzzzz_yzzz[k] = -g_y_0_yyzzz_yzzz[k] * cd_z[k] + g_y_0_yyzzz_yzzzz[k];

                g_y_0_yyzzzz_zzzz[k] = -g_y_0_yyzzz_zzzz[k] * cd_z[k] + g_y_0_yyzzz_zzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 390);

            auto g_y_0_yzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 391);

            auto g_y_0_yzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 392);

            auto g_y_0_yzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 393);

            auto g_y_0_yzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 394);

            auto g_y_0_yzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 395);

            auto g_y_0_yzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 396);

            auto g_y_0_yzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 397);

            auto g_y_0_yzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 398);

            auto g_y_0_yzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 399);

            auto g_y_0_yzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 400);

            auto g_y_0_yzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 401);

            auto g_y_0_yzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 402);

            auto g_y_0_yzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 403);

            auto g_y_0_yzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 404);

            #pragma omp simd aligned(cd_z, g_y_0_yzzzz_xxxx, g_y_0_yzzzz_xxxxz, g_y_0_yzzzz_xxxy, g_y_0_yzzzz_xxxyz, g_y_0_yzzzz_xxxz, g_y_0_yzzzz_xxxzz, g_y_0_yzzzz_xxyy, g_y_0_yzzzz_xxyyz, g_y_0_yzzzz_xxyz, g_y_0_yzzzz_xxyzz, g_y_0_yzzzz_xxzz, g_y_0_yzzzz_xxzzz, g_y_0_yzzzz_xyyy, g_y_0_yzzzz_xyyyz, g_y_0_yzzzz_xyyz, g_y_0_yzzzz_xyyzz, g_y_0_yzzzz_xyzz, g_y_0_yzzzz_xyzzz, g_y_0_yzzzz_xzzz, g_y_0_yzzzz_xzzzz, g_y_0_yzzzz_yyyy, g_y_0_yzzzz_yyyyz, g_y_0_yzzzz_yyyz, g_y_0_yzzzz_yyyzz, g_y_0_yzzzz_yyzz, g_y_0_yzzzz_yyzzz, g_y_0_yzzzz_yzzz, g_y_0_yzzzz_yzzzz, g_y_0_yzzzz_zzzz, g_y_0_yzzzz_zzzzz, g_y_0_yzzzzz_xxxx, g_y_0_yzzzzz_xxxy, g_y_0_yzzzzz_xxxz, g_y_0_yzzzzz_xxyy, g_y_0_yzzzzz_xxyz, g_y_0_yzzzzz_xxzz, g_y_0_yzzzzz_xyyy, g_y_0_yzzzzz_xyyz, g_y_0_yzzzzz_xyzz, g_y_0_yzzzzz_xzzz, g_y_0_yzzzzz_yyyy, g_y_0_yzzzzz_yyyz, g_y_0_yzzzzz_yyzz, g_y_0_yzzzzz_yzzz, g_y_0_yzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_xxxx[k] = -g_y_0_yzzzz_xxxx[k] * cd_z[k] + g_y_0_yzzzz_xxxxz[k];

                g_y_0_yzzzzz_xxxy[k] = -g_y_0_yzzzz_xxxy[k] * cd_z[k] + g_y_0_yzzzz_xxxyz[k];

                g_y_0_yzzzzz_xxxz[k] = -g_y_0_yzzzz_xxxz[k] * cd_z[k] + g_y_0_yzzzz_xxxzz[k];

                g_y_0_yzzzzz_xxyy[k] = -g_y_0_yzzzz_xxyy[k] * cd_z[k] + g_y_0_yzzzz_xxyyz[k];

                g_y_0_yzzzzz_xxyz[k] = -g_y_0_yzzzz_xxyz[k] * cd_z[k] + g_y_0_yzzzz_xxyzz[k];

                g_y_0_yzzzzz_xxzz[k] = -g_y_0_yzzzz_xxzz[k] * cd_z[k] + g_y_0_yzzzz_xxzzz[k];

                g_y_0_yzzzzz_xyyy[k] = -g_y_0_yzzzz_xyyy[k] * cd_z[k] + g_y_0_yzzzz_xyyyz[k];

                g_y_0_yzzzzz_xyyz[k] = -g_y_0_yzzzz_xyyz[k] * cd_z[k] + g_y_0_yzzzz_xyyzz[k];

                g_y_0_yzzzzz_xyzz[k] = -g_y_0_yzzzz_xyzz[k] * cd_z[k] + g_y_0_yzzzz_xyzzz[k];

                g_y_0_yzzzzz_xzzz[k] = -g_y_0_yzzzz_xzzz[k] * cd_z[k] + g_y_0_yzzzz_xzzzz[k];

                g_y_0_yzzzzz_yyyy[k] = -g_y_0_yzzzz_yyyy[k] * cd_z[k] + g_y_0_yzzzz_yyyyz[k];

                g_y_0_yzzzzz_yyyz[k] = -g_y_0_yzzzz_yyyz[k] * cd_z[k] + g_y_0_yzzzz_yyyzz[k];

                g_y_0_yzzzzz_yyzz[k] = -g_y_0_yzzzz_yyzz[k] * cd_z[k] + g_y_0_yzzzz_yyzzz[k];

                g_y_0_yzzzzz_yzzz[k] = -g_y_0_yzzzz_yzzz[k] * cd_z[k] + g_y_0_yzzzz_yzzzz[k];

                g_y_0_yzzzzz_zzzz[k] = -g_y_0_yzzzz_zzzz[k] * cd_z[k] + g_y_0_yzzzz_zzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 405);

            auto g_y_0_zzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 406);

            auto g_y_0_zzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 407);

            auto g_y_0_zzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 408);

            auto g_y_0_zzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 409);

            auto g_y_0_zzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 410);

            auto g_y_0_zzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 411);

            auto g_y_0_zzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 412);

            auto g_y_0_zzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 413);

            auto g_y_0_zzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 414);

            auto g_y_0_zzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 415);

            auto g_y_0_zzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 416);

            auto g_y_0_zzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 417);

            auto g_y_0_zzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 418);

            auto g_y_0_zzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 420 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_z, g_y_0_zzzzz_xxxx, g_y_0_zzzzz_xxxxz, g_y_0_zzzzz_xxxy, g_y_0_zzzzz_xxxyz, g_y_0_zzzzz_xxxz, g_y_0_zzzzz_xxxzz, g_y_0_zzzzz_xxyy, g_y_0_zzzzz_xxyyz, g_y_0_zzzzz_xxyz, g_y_0_zzzzz_xxyzz, g_y_0_zzzzz_xxzz, g_y_0_zzzzz_xxzzz, g_y_0_zzzzz_xyyy, g_y_0_zzzzz_xyyyz, g_y_0_zzzzz_xyyz, g_y_0_zzzzz_xyyzz, g_y_0_zzzzz_xyzz, g_y_0_zzzzz_xyzzz, g_y_0_zzzzz_xzzz, g_y_0_zzzzz_xzzzz, g_y_0_zzzzz_yyyy, g_y_0_zzzzz_yyyyz, g_y_0_zzzzz_yyyz, g_y_0_zzzzz_yyyzz, g_y_0_zzzzz_yyzz, g_y_0_zzzzz_yyzzz, g_y_0_zzzzz_yzzz, g_y_0_zzzzz_yzzzz, g_y_0_zzzzz_zzzz, g_y_0_zzzzz_zzzzz, g_y_0_zzzzzz_xxxx, g_y_0_zzzzzz_xxxy, g_y_0_zzzzzz_xxxz, g_y_0_zzzzzz_xxyy, g_y_0_zzzzzz_xxyz, g_y_0_zzzzzz_xxzz, g_y_0_zzzzzz_xyyy, g_y_0_zzzzzz_xyyz, g_y_0_zzzzzz_xyzz, g_y_0_zzzzzz_xzzz, g_y_0_zzzzzz_yyyy, g_y_0_zzzzzz_yyyz, g_y_0_zzzzzz_yyzz, g_y_0_zzzzzz_yzzz, g_y_0_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_xxxx[k] = -g_y_0_zzzzz_xxxx[k] * cd_z[k] + g_y_0_zzzzz_xxxxz[k];

                g_y_0_zzzzzz_xxxy[k] = -g_y_0_zzzzz_xxxy[k] * cd_z[k] + g_y_0_zzzzz_xxxyz[k];

                g_y_0_zzzzzz_xxxz[k] = -g_y_0_zzzzz_xxxz[k] * cd_z[k] + g_y_0_zzzzz_xxxzz[k];

                g_y_0_zzzzzz_xxyy[k] = -g_y_0_zzzzz_xxyy[k] * cd_z[k] + g_y_0_zzzzz_xxyyz[k];

                g_y_0_zzzzzz_xxyz[k] = -g_y_0_zzzzz_xxyz[k] * cd_z[k] + g_y_0_zzzzz_xxyzz[k];

                g_y_0_zzzzzz_xxzz[k] = -g_y_0_zzzzz_xxzz[k] * cd_z[k] + g_y_0_zzzzz_xxzzz[k];

                g_y_0_zzzzzz_xyyy[k] = -g_y_0_zzzzz_xyyy[k] * cd_z[k] + g_y_0_zzzzz_xyyyz[k];

                g_y_0_zzzzzz_xyyz[k] = -g_y_0_zzzzz_xyyz[k] * cd_z[k] + g_y_0_zzzzz_xyyzz[k];

                g_y_0_zzzzzz_xyzz[k] = -g_y_0_zzzzz_xyzz[k] * cd_z[k] + g_y_0_zzzzz_xyzzz[k];

                g_y_0_zzzzzz_xzzz[k] = -g_y_0_zzzzz_xzzz[k] * cd_z[k] + g_y_0_zzzzz_xzzzz[k];

                g_y_0_zzzzzz_yyyy[k] = -g_y_0_zzzzz_yyyy[k] * cd_z[k] + g_y_0_zzzzz_yyyyz[k];

                g_y_0_zzzzzz_yyyz[k] = -g_y_0_zzzzz_yyyz[k] * cd_z[k] + g_y_0_zzzzz_yyyzz[k];

                g_y_0_zzzzzz_yyzz[k] = -g_y_0_zzzzz_yyzz[k] * cd_z[k] + g_y_0_zzzzz_yyzzz[k];

                g_y_0_zzzzzz_yzzz[k] = -g_y_0_zzzzz_yzzz[k] * cd_z[k] + g_y_0_zzzzz_yzzzz[k];

                g_y_0_zzzzzz_zzzz[k] = -g_y_0_zzzzz_zzzz[k] * cd_z[k] + g_y_0_zzzzz_zzzzz[k];
            }
            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 0);

            auto g_z_0_xxxxxx_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 1);

            auto g_z_0_xxxxxx_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 2);

            auto g_z_0_xxxxxx_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 3);

            auto g_z_0_xxxxxx_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 4);

            auto g_z_0_xxxxxx_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 5);

            auto g_z_0_xxxxxx_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 6);

            auto g_z_0_xxxxxx_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 7);

            auto g_z_0_xxxxxx_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 8);

            auto g_z_0_xxxxxx_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 9);

            auto g_z_0_xxxxxx_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 10);

            auto g_z_0_xxxxxx_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 11);

            auto g_z_0_xxxxxx_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 12);

            auto g_z_0_xxxxxx_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 13);

            auto g_z_0_xxxxxx_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxx_xxxx, g_z_0_xxxxx_xxxxx, g_z_0_xxxxx_xxxxy, g_z_0_xxxxx_xxxxz, g_z_0_xxxxx_xxxy, g_z_0_xxxxx_xxxyy, g_z_0_xxxxx_xxxyz, g_z_0_xxxxx_xxxz, g_z_0_xxxxx_xxxzz, g_z_0_xxxxx_xxyy, g_z_0_xxxxx_xxyyy, g_z_0_xxxxx_xxyyz, g_z_0_xxxxx_xxyz, g_z_0_xxxxx_xxyzz, g_z_0_xxxxx_xxzz, g_z_0_xxxxx_xxzzz, g_z_0_xxxxx_xyyy, g_z_0_xxxxx_xyyyy, g_z_0_xxxxx_xyyyz, g_z_0_xxxxx_xyyz, g_z_0_xxxxx_xyyzz, g_z_0_xxxxx_xyzz, g_z_0_xxxxx_xyzzz, g_z_0_xxxxx_xzzz, g_z_0_xxxxx_xzzzz, g_z_0_xxxxx_yyyy, g_z_0_xxxxx_yyyz, g_z_0_xxxxx_yyzz, g_z_0_xxxxx_yzzz, g_z_0_xxxxx_zzzz, g_z_0_xxxxxx_xxxx, g_z_0_xxxxxx_xxxy, g_z_0_xxxxxx_xxxz, g_z_0_xxxxxx_xxyy, g_z_0_xxxxxx_xxyz, g_z_0_xxxxxx_xxzz, g_z_0_xxxxxx_xyyy, g_z_0_xxxxxx_xyyz, g_z_0_xxxxxx_xyzz, g_z_0_xxxxxx_xzzz, g_z_0_xxxxxx_yyyy, g_z_0_xxxxxx_yyyz, g_z_0_xxxxxx_yyzz, g_z_0_xxxxxx_yzzz, g_z_0_xxxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_xxxx[k] = -g_z_0_xxxxx_xxxx[k] * cd_x[k] + g_z_0_xxxxx_xxxxx[k];

                g_z_0_xxxxxx_xxxy[k] = -g_z_0_xxxxx_xxxy[k] * cd_x[k] + g_z_0_xxxxx_xxxxy[k];

                g_z_0_xxxxxx_xxxz[k] = -g_z_0_xxxxx_xxxz[k] * cd_x[k] + g_z_0_xxxxx_xxxxz[k];

                g_z_0_xxxxxx_xxyy[k] = -g_z_0_xxxxx_xxyy[k] * cd_x[k] + g_z_0_xxxxx_xxxyy[k];

                g_z_0_xxxxxx_xxyz[k] = -g_z_0_xxxxx_xxyz[k] * cd_x[k] + g_z_0_xxxxx_xxxyz[k];

                g_z_0_xxxxxx_xxzz[k] = -g_z_0_xxxxx_xxzz[k] * cd_x[k] + g_z_0_xxxxx_xxxzz[k];

                g_z_0_xxxxxx_xyyy[k] = -g_z_0_xxxxx_xyyy[k] * cd_x[k] + g_z_0_xxxxx_xxyyy[k];

                g_z_0_xxxxxx_xyyz[k] = -g_z_0_xxxxx_xyyz[k] * cd_x[k] + g_z_0_xxxxx_xxyyz[k];

                g_z_0_xxxxxx_xyzz[k] = -g_z_0_xxxxx_xyzz[k] * cd_x[k] + g_z_0_xxxxx_xxyzz[k];

                g_z_0_xxxxxx_xzzz[k] = -g_z_0_xxxxx_xzzz[k] * cd_x[k] + g_z_0_xxxxx_xxzzz[k];

                g_z_0_xxxxxx_yyyy[k] = -g_z_0_xxxxx_yyyy[k] * cd_x[k] + g_z_0_xxxxx_xyyyy[k];

                g_z_0_xxxxxx_yyyz[k] = -g_z_0_xxxxx_yyyz[k] * cd_x[k] + g_z_0_xxxxx_xyyyz[k];

                g_z_0_xxxxxx_yyzz[k] = -g_z_0_xxxxx_yyzz[k] * cd_x[k] + g_z_0_xxxxx_xyyzz[k];

                g_z_0_xxxxxx_yzzz[k] = -g_z_0_xxxxx_yzzz[k] * cd_x[k] + g_z_0_xxxxx_xyzzz[k];

                g_z_0_xxxxxx_zzzz[k] = -g_z_0_xxxxx_zzzz[k] * cd_x[k] + g_z_0_xxxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 15);

            auto g_z_0_xxxxxy_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 16);

            auto g_z_0_xxxxxy_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 17);

            auto g_z_0_xxxxxy_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 18);

            auto g_z_0_xxxxxy_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 19);

            auto g_z_0_xxxxxy_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 20);

            auto g_z_0_xxxxxy_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 21);

            auto g_z_0_xxxxxy_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 22);

            auto g_z_0_xxxxxy_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 23);

            auto g_z_0_xxxxxy_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 24);

            auto g_z_0_xxxxxy_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 25);

            auto g_z_0_xxxxxy_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 26);

            auto g_z_0_xxxxxy_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 27);

            auto g_z_0_xxxxxy_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 28);

            auto g_z_0_xxxxxy_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxy_xxxx, g_z_0_xxxxxy_xxxy, g_z_0_xxxxxy_xxxz, g_z_0_xxxxxy_xxyy, g_z_0_xxxxxy_xxyz, g_z_0_xxxxxy_xxzz, g_z_0_xxxxxy_xyyy, g_z_0_xxxxxy_xyyz, g_z_0_xxxxxy_xyzz, g_z_0_xxxxxy_xzzz, g_z_0_xxxxxy_yyyy, g_z_0_xxxxxy_yyyz, g_z_0_xxxxxy_yyzz, g_z_0_xxxxxy_yzzz, g_z_0_xxxxxy_zzzz, g_z_0_xxxxy_xxxx, g_z_0_xxxxy_xxxxx, g_z_0_xxxxy_xxxxy, g_z_0_xxxxy_xxxxz, g_z_0_xxxxy_xxxy, g_z_0_xxxxy_xxxyy, g_z_0_xxxxy_xxxyz, g_z_0_xxxxy_xxxz, g_z_0_xxxxy_xxxzz, g_z_0_xxxxy_xxyy, g_z_0_xxxxy_xxyyy, g_z_0_xxxxy_xxyyz, g_z_0_xxxxy_xxyz, g_z_0_xxxxy_xxyzz, g_z_0_xxxxy_xxzz, g_z_0_xxxxy_xxzzz, g_z_0_xxxxy_xyyy, g_z_0_xxxxy_xyyyy, g_z_0_xxxxy_xyyyz, g_z_0_xxxxy_xyyz, g_z_0_xxxxy_xyyzz, g_z_0_xxxxy_xyzz, g_z_0_xxxxy_xyzzz, g_z_0_xxxxy_xzzz, g_z_0_xxxxy_xzzzz, g_z_0_xxxxy_yyyy, g_z_0_xxxxy_yyyz, g_z_0_xxxxy_yyzz, g_z_0_xxxxy_yzzz, g_z_0_xxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_xxxx[k] = -g_z_0_xxxxy_xxxx[k] * cd_x[k] + g_z_0_xxxxy_xxxxx[k];

                g_z_0_xxxxxy_xxxy[k] = -g_z_0_xxxxy_xxxy[k] * cd_x[k] + g_z_0_xxxxy_xxxxy[k];

                g_z_0_xxxxxy_xxxz[k] = -g_z_0_xxxxy_xxxz[k] * cd_x[k] + g_z_0_xxxxy_xxxxz[k];

                g_z_0_xxxxxy_xxyy[k] = -g_z_0_xxxxy_xxyy[k] * cd_x[k] + g_z_0_xxxxy_xxxyy[k];

                g_z_0_xxxxxy_xxyz[k] = -g_z_0_xxxxy_xxyz[k] * cd_x[k] + g_z_0_xxxxy_xxxyz[k];

                g_z_0_xxxxxy_xxzz[k] = -g_z_0_xxxxy_xxzz[k] * cd_x[k] + g_z_0_xxxxy_xxxzz[k];

                g_z_0_xxxxxy_xyyy[k] = -g_z_0_xxxxy_xyyy[k] * cd_x[k] + g_z_0_xxxxy_xxyyy[k];

                g_z_0_xxxxxy_xyyz[k] = -g_z_0_xxxxy_xyyz[k] * cd_x[k] + g_z_0_xxxxy_xxyyz[k];

                g_z_0_xxxxxy_xyzz[k] = -g_z_0_xxxxy_xyzz[k] * cd_x[k] + g_z_0_xxxxy_xxyzz[k];

                g_z_0_xxxxxy_xzzz[k] = -g_z_0_xxxxy_xzzz[k] * cd_x[k] + g_z_0_xxxxy_xxzzz[k];

                g_z_0_xxxxxy_yyyy[k] = -g_z_0_xxxxy_yyyy[k] * cd_x[k] + g_z_0_xxxxy_xyyyy[k];

                g_z_0_xxxxxy_yyyz[k] = -g_z_0_xxxxy_yyyz[k] * cd_x[k] + g_z_0_xxxxy_xyyyz[k];

                g_z_0_xxxxxy_yyzz[k] = -g_z_0_xxxxy_yyzz[k] * cd_x[k] + g_z_0_xxxxy_xyyzz[k];

                g_z_0_xxxxxy_yzzz[k] = -g_z_0_xxxxy_yzzz[k] * cd_x[k] + g_z_0_xxxxy_xyzzz[k];

                g_z_0_xxxxxy_zzzz[k] = -g_z_0_xxxxy_zzzz[k] * cd_x[k] + g_z_0_xxxxy_xzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 30);

            auto g_z_0_xxxxxz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 31);

            auto g_z_0_xxxxxz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 32);

            auto g_z_0_xxxxxz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 33);

            auto g_z_0_xxxxxz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 34);

            auto g_z_0_xxxxxz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 35);

            auto g_z_0_xxxxxz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 36);

            auto g_z_0_xxxxxz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 37);

            auto g_z_0_xxxxxz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 38);

            auto g_z_0_xxxxxz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 39);

            auto g_z_0_xxxxxz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 40);

            auto g_z_0_xxxxxz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 41);

            auto g_z_0_xxxxxz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 42);

            auto g_z_0_xxxxxz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 43);

            auto g_z_0_xxxxxz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxz_xxxx, g_z_0_xxxxxz_xxxy, g_z_0_xxxxxz_xxxz, g_z_0_xxxxxz_xxyy, g_z_0_xxxxxz_xxyz, g_z_0_xxxxxz_xxzz, g_z_0_xxxxxz_xyyy, g_z_0_xxxxxz_xyyz, g_z_0_xxxxxz_xyzz, g_z_0_xxxxxz_xzzz, g_z_0_xxxxxz_yyyy, g_z_0_xxxxxz_yyyz, g_z_0_xxxxxz_yyzz, g_z_0_xxxxxz_yzzz, g_z_0_xxxxxz_zzzz, g_z_0_xxxxz_xxxx, g_z_0_xxxxz_xxxxx, g_z_0_xxxxz_xxxxy, g_z_0_xxxxz_xxxxz, g_z_0_xxxxz_xxxy, g_z_0_xxxxz_xxxyy, g_z_0_xxxxz_xxxyz, g_z_0_xxxxz_xxxz, g_z_0_xxxxz_xxxzz, g_z_0_xxxxz_xxyy, g_z_0_xxxxz_xxyyy, g_z_0_xxxxz_xxyyz, g_z_0_xxxxz_xxyz, g_z_0_xxxxz_xxyzz, g_z_0_xxxxz_xxzz, g_z_0_xxxxz_xxzzz, g_z_0_xxxxz_xyyy, g_z_0_xxxxz_xyyyy, g_z_0_xxxxz_xyyyz, g_z_0_xxxxz_xyyz, g_z_0_xxxxz_xyyzz, g_z_0_xxxxz_xyzz, g_z_0_xxxxz_xyzzz, g_z_0_xxxxz_xzzz, g_z_0_xxxxz_xzzzz, g_z_0_xxxxz_yyyy, g_z_0_xxxxz_yyyz, g_z_0_xxxxz_yyzz, g_z_0_xxxxz_yzzz, g_z_0_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_xxxx[k] = -g_z_0_xxxxz_xxxx[k] * cd_x[k] + g_z_0_xxxxz_xxxxx[k];

                g_z_0_xxxxxz_xxxy[k] = -g_z_0_xxxxz_xxxy[k] * cd_x[k] + g_z_0_xxxxz_xxxxy[k];

                g_z_0_xxxxxz_xxxz[k] = -g_z_0_xxxxz_xxxz[k] * cd_x[k] + g_z_0_xxxxz_xxxxz[k];

                g_z_0_xxxxxz_xxyy[k] = -g_z_0_xxxxz_xxyy[k] * cd_x[k] + g_z_0_xxxxz_xxxyy[k];

                g_z_0_xxxxxz_xxyz[k] = -g_z_0_xxxxz_xxyz[k] * cd_x[k] + g_z_0_xxxxz_xxxyz[k];

                g_z_0_xxxxxz_xxzz[k] = -g_z_0_xxxxz_xxzz[k] * cd_x[k] + g_z_0_xxxxz_xxxzz[k];

                g_z_0_xxxxxz_xyyy[k] = -g_z_0_xxxxz_xyyy[k] * cd_x[k] + g_z_0_xxxxz_xxyyy[k];

                g_z_0_xxxxxz_xyyz[k] = -g_z_0_xxxxz_xyyz[k] * cd_x[k] + g_z_0_xxxxz_xxyyz[k];

                g_z_0_xxxxxz_xyzz[k] = -g_z_0_xxxxz_xyzz[k] * cd_x[k] + g_z_0_xxxxz_xxyzz[k];

                g_z_0_xxxxxz_xzzz[k] = -g_z_0_xxxxz_xzzz[k] * cd_x[k] + g_z_0_xxxxz_xxzzz[k];

                g_z_0_xxxxxz_yyyy[k] = -g_z_0_xxxxz_yyyy[k] * cd_x[k] + g_z_0_xxxxz_xyyyy[k];

                g_z_0_xxxxxz_yyyz[k] = -g_z_0_xxxxz_yyyz[k] * cd_x[k] + g_z_0_xxxxz_xyyyz[k];

                g_z_0_xxxxxz_yyzz[k] = -g_z_0_xxxxz_yyzz[k] * cd_x[k] + g_z_0_xxxxz_xyyzz[k];

                g_z_0_xxxxxz_yzzz[k] = -g_z_0_xxxxz_yzzz[k] * cd_x[k] + g_z_0_xxxxz_xyzzz[k];

                g_z_0_xxxxxz_zzzz[k] = -g_z_0_xxxxz_zzzz[k] * cd_x[k] + g_z_0_xxxxz_xzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 45);

            auto g_z_0_xxxxyy_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 46);

            auto g_z_0_xxxxyy_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 47);

            auto g_z_0_xxxxyy_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 48);

            auto g_z_0_xxxxyy_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 49);

            auto g_z_0_xxxxyy_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 50);

            auto g_z_0_xxxxyy_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 51);

            auto g_z_0_xxxxyy_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 52);

            auto g_z_0_xxxxyy_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 53);

            auto g_z_0_xxxxyy_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 54);

            auto g_z_0_xxxxyy_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 55);

            auto g_z_0_xxxxyy_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 56);

            auto g_z_0_xxxxyy_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 57);

            auto g_z_0_xxxxyy_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 58);

            auto g_z_0_xxxxyy_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyy_xxxx, g_z_0_xxxxyy_xxxy, g_z_0_xxxxyy_xxxz, g_z_0_xxxxyy_xxyy, g_z_0_xxxxyy_xxyz, g_z_0_xxxxyy_xxzz, g_z_0_xxxxyy_xyyy, g_z_0_xxxxyy_xyyz, g_z_0_xxxxyy_xyzz, g_z_0_xxxxyy_xzzz, g_z_0_xxxxyy_yyyy, g_z_0_xxxxyy_yyyz, g_z_0_xxxxyy_yyzz, g_z_0_xxxxyy_yzzz, g_z_0_xxxxyy_zzzz, g_z_0_xxxyy_xxxx, g_z_0_xxxyy_xxxxx, g_z_0_xxxyy_xxxxy, g_z_0_xxxyy_xxxxz, g_z_0_xxxyy_xxxy, g_z_0_xxxyy_xxxyy, g_z_0_xxxyy_xxxyz, g_z_0_xxxyy_xxxz, g_z_0_xxxyy_xxxzz, g_z_0_xxxyy_xxyy, g_z_0_xxxyy_xxyyy, g_z_0_xxxyy_xxyyz, g_z_0_xxxyy_xxyz, g_z_0_xxxyy_xxyzz, g_z_0_xxxyy_xxzz, g_z_0_xxxyy_xxzzz, g_z_0_xxxyy_xyyy, g_z_0_xxxyy_xyyyy, g_z_0_xxxyy_xyyyz, g_z_0_xxxyy_xyyz, g_z_0_xxxyy_xyyzz, g_z_0_xxxyy_xyzz, g_z_0_xxxyy_xyzzz, g_z_0_xxxyy_xzzz, g_z_0_xxxyy_xzzzz, g_z_0_xxxyy_yyyy, g_z_0_xxxyy_yyyz, g_z_0_xxxyy_yyzz, g_z_0_xxxyy_yzzz, g_z_0_xxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_xxxx[k] = -g_z_0_xxxyy_xxxx[k] * cd_x[k] + g_z_0_xxxyy_xxxxx[k];

                g_z_0_xxxxyy_xxxy[k] = -g_z_0_xxxyy_xxxy[k] * cd_x[k] + g_z_0_xxxyy_xxxxy[k];

                g_z_0_xxxxyy_xxxz[k] = -g_z_0_xxxyy_xxxz[k] * cd_x[k] + g_z_0_xxxyy_xxxxz[k];

                g_z_0_xxxxyy_xxyy[k] = -g_z_0_xxxyy_xxyy[k] * cd_x[k] + g_z_0_xxxyy_xxxyy[k];

                g_z_0_xxxxyy_xxyz[k] = -g_z_0_xxxyy_xxyz[k] * cd_x[k] + g_z_0_xxxyy_xxxyz[k];

                g_z_0_xxxxyy_xxzz[k] = -g_z_0_xxxyy_xxzz[k] * cd_x[k] + g_z_0_xxxyy_xxxzz[k];

                g_z_0_xxxxyy_xyyy[k] = -g_z_0_xxxyy_xyyy[k] * cd_x[k] + g_z_0_xxxyy_xxyyy[k];

                g_z_0_xxxxyy_xyyz[k] = -g_z_0_xxxyy_xyyz[k] * cd_x[k] + g_z_0_xxxyy_xxyyz[k];

                g_z_0_xxxxyy_xyzz[k] = -g_z_0_xxxyy_xyzz[k] * cd_x[k] + g_z_0_xxxyy_xxyzz[k];

                g_z_0_xxxxyy_xzzz[k] = -g_z_0_xxxyy_xzzz[k] * cd_x[k] + g_z_0_xxxyy_xxzzz[k];

                g_z_0_xxxxyy_yyyy[k] = -g_z_0_xxxyy_yyyy[k] * cd_x[k] + g_z_0_xxxyy_xyyyy[k];

                g_z_0_xxxxyy_yyyz[k] = -g_z_0_xxxyy_yyyz[k] * cd_x[k] + g_z_0_xxxyy_xyyyz[k];

                g_z_0_xxxxyy_yyzz[k] = -g_z_0_xxxyy_yyzz[k] * cd_x[k] + g_z_0_xxxyy_xyyzz[k];

                g_z_0_xxxxyy_yzzz[k] = -g_z_0_xxxyy_yzzz[k] * cd_x[k] + g_z_0_xxxyy_xyzzz[k];

                g_z_0_xxxxyy_zzzz[k] = -g_z_0_xxxyy_zzzz[k] * cd_x[k] + g_z_0_xxxyy_xzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 60);

            auto g_z_0_xxxxyz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 61);

            auto g_z_0_xxxxyz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 62);

            auto g_z_0_xxxxyz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 63);

            auto g_z_0_xxxxyz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 64);

            auto g_z_0_xxxxyz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 65);

            auto g_z_0_xxxxyz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 66);

            auto g_z_0_xxxxyz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 67);

            auto g_z_0_xxxxyz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 68);

            auto g_z_0_xxxxyz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 69);

            auto g_z_0_xxxxyz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 70);

            auto g_z_0_xxxxyz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 71);

            auto g_z_0_xxxxyz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 72);

            auto g_z_0_xxxxyz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 73);

            auto g_z_0_xxxxyz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 74);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyz_xxxx, g_z_0_xxxxyz_xxxy, g_z_0_xxxxyz_xxxz, g_z_0_xxxxyz_xxyy, g_z_0_xxxxyz_xxyz, g_z_0_xxxxyz_xxzz, g_z_0_xxxxyz_xyyy, g_z_0_xxxxyz_xyyz, g_z_0_xxxxyz_xyzz, g_z_0_xxxxyz_xzzz, g_z_0_xxxxyz_yyyy, g_z_0_xxxxyz_yyyz, g_z_0_xxxxyz_yyzz, g_z_0_xxxxyz_yzzz, g_z_0_xxxxyz_zzzz, g_z_0_xxxyz_xxxx, g_z_0_xxxyz_xxxxx, g_z_0_xxxyz_xxxxy, g_z_0_xxxyz_xxxxz, g_z_0_xxxyz_xxxy, g_z_0_xxxyz_xxxyy, g_z_0_xxxyz_xxxyz, g_z_0_xxxyz_xxxz, g_z_0_xxxyz_xxxzz, g_z_0_xxxyz_xxyy, g_z_0_xxxyz_xxyyy, g_z_0_xxxyz_xxyyz, g_z_0_xxxyz_xxyz, g_z_0_xxxyz_xxyzz, g_z_0_xxxyz_xxzz, g_z_0_xxxyz_xxzzz, g_z_0_xxxyz_xyyy, g_z_0_xxxyz_xyyyy, g_z_0_xxxyz_xyyyz, g_z_0_xxxyz_xyyz, g_z_0_xxxyz_xyyzz, g_z_0_xxxyz_xyzz, g_z_0_xxxyz_xyzzz, g_z_0_xxxyz_xzzz, g_z_0_xxxyz_xzzzz, g_z_0_xxxyz_yyyy, g_z_0_xxxyz_yyyz, g_z_0_xxxyz_yyzz, g_z_0_xxxyz_yzzz, g_z_0_xxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_xxxx[k] = -g_z_0_xxxyz_xxxx[k] * cd_x[k] + g_z_0_xxxyz_xxxxx[k];

                g_z_0_xxxxyz_xxxy[k] = -g_z_0_xxxyz_xxxy[k] * cd_x[k] + g_z_0_xxxyz_xxxxy[k];

                g_z_0_xxxxyz_xxxz[k] = -g_z_0_xxxyz_xxxz[k] * cd_x[k] + g_z_0_xxxyz_xxxxz[k];

                g_z_0_xxxxyz_xxyy[k] = -g_z_0_xxxyz_xxyy[k] * cd_x[k] + g_z_0_xxxyz_xxxyy[k];

                g_z_0_xxxxyz_xxyz[k] = -g_z_0_xxxyz_xxyz[k] * cd_x[k] + g_z_0_xxxyz_xxxyz[k];

                g_z_0_xxxxyz_xxzz[k] = -g_z_0_xxxyz_xxzz[k] * cd_x[k] + g_z_0_xxxyz_xxxzz[k];

                g_z_0_xxxxyz_xyyy[k] = -g_z_0_xxxyz_xyyy[k] * cd_x[k] + g_z_0_xxxyz_xxyyy[k];

                g_z_0_xxxxyz_xyyz[k] = -g_z_0_xxxyz_xyyz[k] * cd_x[k] + g_z_0_xxxyz_xxyyz[k];

                g_z_0_xxxxyz_xyzz[k] = -g_z_0_xxxyz_xyzz[k] * cd_x[k] + g_z_0_xxxyz_xxyzz[k];

                g_z_0_xxxxyz_xzzz[k] = -g_z_0_xxxyz_xzzz[k] * cd_x[k] + g_z_0_xxxyz_xxzzz[k];

                g_z_0_xxxxyz_yyyy[k] = -g_z_0_xxxyz_yyyy[k] * cd_x[k] + g_z_0_xxxyz_xyyyy[k];

                g_z_0_xxxxyz_yyyz[k] = -g_z_0_xxxyz_yyyz[k] * cd_x[k] + g_z_0_xxxyz_xyyyz[k];

                g_z_0_xxxxyz_yyzz[k] = -g_z_0_xxxyz_yyzz[k] * cd_x[k] + g_z_0_xxxyz_xyyzz[k];

                g_z_0_xxxxyz_yzzz[k] = -g_z_0_xxxyz_yzzz[k] * cd_x[k] + g_z_0_xxxyz_xyzzz[k];

                g_z_0_xxxxyz_zzzz[k] = -g_z_0_xxxyz_zzzz[k] * cd_x[k] + g_z_0_xxxyz_xzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 75);

            auto g_z_0_xxxxzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 76);

            auto g_z_0_xxxxzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 77);

            auto g_z_0_xxxxzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 78);

            auto g_z_0_xxxxzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 79);

            auto g_z_0_xxxxzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 80);

            auto g_z_0_xxxxzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 81);

            auto g_z_0_xxxxzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 82);

            auto g_z_0_xxxxzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 83);

            auto g_z_0_xxxxzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 84);

            auto g_z_0_xxxxzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 85);

            auto g_z_0_xxxxzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 86);

            auto g_z_0_xxxxzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 87);

            auto g_z_0_xxxxzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 88);

            auto g_z_0_xxxxzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxzz_xxxx, g_z_0_xxxxzz_xxxy, g_z_0_xxxxzz_xxxz, g_z_0_xxxxzz_xxyy, g_z_0_xxxxzz_xxyz, g_z_0_xxxxzz_xxzz, g_z_0_xxxxzz_xyyy, g_z_0_xxxxzz_xyyz, g_z_0_xxxxzz_xyzz, g_z_0_xxxxzz_xzzz, g_z_0_xxxxzz_yyyy, g_z_0_xxxxzz_yyyz, g_z_0_xxxxzz_yyzz, g_z_0_xxxxzz_yzzz, g_z_0_xxxxzz_zzzz, g_z_0_xxxzz_xxxx, g_z_0_xxxzz_xxxxx, g_z_0_xxxzz_xxxxy, g_z_0_xxxzz_xxxxz, g_z_0_xxxzz_xxxy, g_z_0_xxxzz_xxxyy, g_z_0_xxxzz_xxxyz, g_z_0_xxxzz_xxxz, g_z_0_xxxzz_xxxzz, g_z_0_xxxzz_xxyy, g_z_0_xxxzz_xxyyy, g_z_0_xxxzz_xxyyz, g_z_0_xxxzz_xxyz, g_z_0_xxxzz_xxyzz, g_z_0_xxxzz_xxzz, g_z_0_xxxzz_xxzzz, g_z_0_xxxzz_xyyy, g_z_0_xxxzz_xyyyy, g_z_0_xxxzz_xyyyz, g_z_0_xxxzz_xyyz, g_z_0_xxxzz_xyyzz, g_z_0_xxxzz_xyzz, g_z_0_xxxzz_xyzzz, g_z_0_xxxzz_xzzz, g_z_0_xxxzz_xzzzz, g_z_0_xxxzz_yyyy, g_z_0_xxxzz_yyyz, g_z_0_xxxzz_yyzz, g_z_0_xxxzz_yzzz, g_z_0_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_xxxx[k] = -g_z_0_xxxzz_xxxx[k] * cd_x[k] + g_z_0_xxxzz_xxxxx[k];

                g_z_0_xxxxzz_xxxy[k] = -g_z_0_xxxzz_xxxy[k] * cd_x[k] + g_z_0_xxxzz_xxxxy[k];

                g_z_0_xxxxzz_xxxz[k] = -g_z_0_xxxzz_xxxz[k] * cd_x[k] + g_z_0_xxxzz_xxxxz[k];

                g_z_0_xxxxzz_xxyy[k] = -g_z_0_xxxzz_xxyy[k] * cd_x[k] + g_z_0_xxxzz_xxxyy[k];

                g_z_0_xxxxzz_xxyz[k] = -g_z_0_xxxzz_xxyz[k] * cd_x[k] + g_z_0_xxxzz_xxxyz[k];

                g_z_0_xxxxzz_xxzz[k] = -g_z_0_xxxzz_xxzz[k] * cd_x[k] + g_z_0_xxxzz_xxxzz[k];

                g_z_0_xxxxzz_xyyy[k] = -g_z_0_xxxzz_xyyy[k] * cd_x[k] + g_z_0_xxxzz_xxyyy[k];

                g_z_0_xxxxzz_xyyz[k] = -g_z_0_xxxzz_xyyz[k] * cd_x[k] + g_z_0_xxxzz_xxyyz[k];

                g_z_0_xxxxzz_xyzz[k] = -g_z_0_xxxzz_xyzz[k] * cd_x[k] + g_z_0_xxxzz_xxyzz[k];

                g_z_0_xxxxzz_xzzz[k] = -g_z_0_xxxzz_xzzz[k] * cd_x[k] + g_z_0_xxxzz_xxzzz[k];

                g_z_0_xxxxzz_yyyy[k] = -g_z_0_xxxzz_yyyy[k] * cd_x[k] + g_z_0_xxxzz_xyyyy[k];

                g_z_0_xxxxzz_yyyz[k] = -g_z_0_xxxzz_yyyz[k] * cd_x[k] + g_z_0_xxxzz_xyyyz[k];

                g_z_0_xxxxzz_yyzz[k] = -g_z_0_xxxzz_yyzz[k] * cd_x[k] + g_z_0_xxxzz_xyyzz[k];

                g_z_0_xxxxzz_yzzz[k] = -g_z_0_xxxzz_yzzz[k] * cd_x[k] + g_z_0_xxxzz_xyzzz[k];

                g_z_0_xxxxzz_zzzz[k] = -g_z_0_xxxzz_zzzz[k] * cd_x[k] + g_z_0_xxxzz_xzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 90);

            auto g_z_0_xxxyyy_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 91);

            auto g_z_0_xxxyyy_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 92);

            auto g_z_0_xxxyyy_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 93);

            auto g_z_0_xxxyyy_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 94);

            auto g_z_0_xxxyyy_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 95);

            auto g_z_0_xxxyyy_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 96);

            auto g_z_0_xxxyyy_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 97);

            auto g_z_0_xxxyyy_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 98);

            auto g_z_0_xxxyyy_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 99);

            auto g_z_0_xxxyyy_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 100);

            auto g_z_0_xxxyyy_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 101);

            auto g_z_0_xxxyyy_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 102);

            auto g_z_0_xxxyyy_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 103);

            auto g_z_0_xxxyyy_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyy_xxxx, g_z_0_xxxyyy_xxxy, g_z_0_xxxyyy_xxxz, g_z_0_xxxyyy_xxyy, g_z_0_xxxyyy_xxyz, g_z_0_xxxyyy_xxzz, g_z_0_xxxyyy_xyyy, g_z_0_xxxyyy_xyyz, g_z_0_xxxyyy_xyzz, g_z_0_xxxyyy_xzzz, g_z_0_xxxyyy_yyyy, g_z_0_xxxyyy_yyyz, g_z_0_xxxyyy_yyzz, g_z_0_xxxyyy_yzzz, g_z_0_xxxyyy_zzzz, g_z_0_xxyyy_xxxx, g_z_0_xxyyy_xxxxx, g_z_0_xxyyy_xxxxy, g_z_0_xxyyy_xxxxz, g_z_0_xxyyy_xxxy, g_z_0_xxyyy_xxxyy, g_z_0_xxyyy_xxxyz, g_z_0_xxyyy_xxxz, g_z_0_xxyyy_xxxzz, g_z_0_xxyyy_xxyy, g_z_0_xxyyy_xxyyy, g_z_0_xxyyy_xxyyz, g_z_0_xxyyy_xxyz, g_z_0_xxyyy_xxyzz, g_z_0_xxyyy_xxzz, g_z_0_xxyyy_xxzzz, g_z_0_xxyyy_xyyy, g_z_0_xxyyy_xyyyy, g_z_0_xxyyy_xyyyz, g_z_0_xxyyy_xyyz, g_z_0_xxyyy_xyyzz, g_z_0_xxyyy_xyzz, g_z_0_xxyyy_xyzzz, g_z_0_xxyyy_xzzz, g_z_0_xxyyy_xzzzz, g_z_0_xxyyy_yyyy, g_z_0_xxyyy_yyyz, g_z_0_xxyyy_yyzz, g_z_0_xxyyy_yzzz, g_z_0_xxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_xxxx[k] = -g_z_0_xxyyy_xxxx[k] * cd_x[k] + g_z_0_xxyyy_xxxxx[k];

                g_z_0_xxxyyy_xxxy[k] = -g_z_0_xxyyy_xxxy[k] * cd_x[k] + g_z_0_xxyyy_xxxxy[k];

                g_z_0_xxxyyy_xxxz[k] = -g_z_0_xxyyy_xxxz[k] * cd_x[k] + g_z_0_xxyyy_xxxxz[k];

                g_z_0_xxxyyy_xxyy[k] = -g_z_0_xxyyy_xxyy[k] * cd_x[k] + g_z_0_xxyyy_xxxyy[k];

                g_z_0_xxxyyy_xxyz[k] = -g_z_0_xxyyy_xxyz[k] * cd_x[k] + g_z_0_xxyyy_xxxyz[k];

                g_z_0_xxxyyy_xxzz[k] = -g_z_0_xxyyy_xxzz[k] * cd_x[k] + g_z_0_xxyyy_xxxzz[k];

                g_z_0_xxxyyy_xyyy[k] = -g_z_0_xxyyy_xyyy[k] * cd_x[k] + g_z_0_xxyyy_xxyyy[k];

                g_z_0_xxxyyy_xyyz[k] = -g_z_0_xxyyy_xyyz[k] * cd_x[k] + g_z_0_xxyyy_xxyyz[k];

                g_z_0_xxxyyy_xyzz[k] = -g_z_0_xxyyy_xyzz[k] * cd_x[k] + g_z_0_xxyyy_xxyzz[k];

                g_z_0_xxxyyy_xzzz[k] = -g_z_0_xxyyy_xzzz[k] * cd_x[k] + g_z_0_xxyyy_xxzzz[k];

                g_z_0_xxxyyy_yyyy[k] = -g_z_0_xxyyy_yyyy[k] * cd_x[k] + g_z_0_xxyyy_xyyyy[k];

                g_z_0_xxxyyy_yyyz[k] = -g_z_0_xxyyy_yyyz[k] * cd_x[k] + g_z_0_xxyyy_xyyyz[k];

                g_z_0_xxxyyy_yyzz[k] = -g_z_0_xxyyy_yyzz[k] * cd_x[k] + g_z_0_xxyyy_xyyzz[k];

                g_z_0_xxxyyy_yzzz[k] = -g_z_0_xxyyy_yzzz[k] * cd_x[k] + g_z_0_xxyyy_xyzzz[k];

                g_z_0_xxxyyy_zzzz[k] = -g_z_0_xxyyy_zzzz[k] * cd_x[k] + g_z_0_xxyyy_xzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 105);

            auto g_z_0_xxxyyz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 106);

            auto g_z_0_xxxyyz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 107);

            auto g_z_0_xxxyyz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 108);

            auto g_z_0_xxxyyz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 109);

            auto g_z_0_xxxyyz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 110);

            auto g_z_0_xxxyyz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 111);

            auto g_z_0_xxxyyz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 112);

            auto g_z_0_xxxyyz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 113);

            auto g_z_0_xxxyyz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 114);

            auto g_z_0_xxxyyz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 115);

            auto g_z_0_xxxyyz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 116);

            auto g_z_0_xxxyyz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 117);

            auto g_z_0_xxxyyz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 118);

            auto g_z_0_xxxyyz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyz_xxxx, g_z_0_xxxyyz_xxxy, g_z_0_xxxyyz_xxxz, g_z_0_xxxyyz_xxyy, g_z_0_xxxyyz_xxyz, g_z_0_xxxyyz_xxzz, g_z_0_xxxyyz_xyyy, g_z_0_xxxyyz_xyyz, g_z_0_xxxyyz_xyzz, g_z_0_xxxyyz_xzzz, g_z_0_xxxyyz_yyyy, g_z_0_xxxyyz_yyyz, g_z_0_xxxyyz_yyzz, g_z_0_xxxyyz_yzzz, g_z_0_xxxyyz_zzzz, g_z_0_xxyyz_xxxx, g_z_0_xxyyz_xxxxx, g_z_0_xxyyz_xxxxy, g_z_0_xxyyz_xxxxz, g_z_0_xxyyz_xxxy, g_z_0_xxyyz_xxxyy, g_z_0_xxyyz_xxxyz, g_z_0_xxyyz_xxxz, g_z_0_xxyyz_xxxzz, g_z_0_xxyyz_xxyy, g_z_0_xxyyz_xxyyy, g_z_0_xxyyz_xxyyz, g_z_0_xxyyz_xxyz, g_z_0_xxyyz_xxyzz, g_z_0_xxyyz_xxzz, g_z_0_xxyyz_xxzzz, g_z_0_xxyyz_xyyy, g_z_0_xxyyz_xyyyy, g_z_0_xxyyz_xyyyz, g_z_0_xxyyz_xyyz, g_z_0_xxyyz_xyyzz, g_z_0_xxyyz_xyzz, g_z_0_xxyyz_xyzzz, g_z_0_xxyyz_xzzz, g_z_0_xxyyz_xzzzz, g_z_0_xxyyz_yyyy, g_z_0_xxyyz_yyyz, g_z_0_xxyyz_yyzz, g_z_0_xxyyz_yzzz, g_z_0_xxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_xxxx[k] = -g_z_0_xxyyz_xxxx[k] * cd_x[k] + g_z_0_xxyyz_xxxxx[k];

                g_z_0_xxxyyz_xxxy[k] = -g_z_0_xxyyz_xxxy[k] * cd_x[k] + g_z_0_xxyyz_xxxxy[k];

                g_z_0_xxxyyz_xxxz[k] = -g_z_0_xxyyz_xxxz[k] * cd_x[k] + g_z_0_xxyyz_xxxxz[k];

                g_z_0_xxxyyz_xxyy[k] = -g_z_0_xxyyz_xxyy[k] * cd_x[k] + g_z_0_xxyyz_xxxyy[k];

                g_z_0_xxxyyz_xxyz[k] = -g_z_0_xxyyz_xxyz[k] * cd_x[k] + g_z_0_xxyyz_xxxyz[k];

                g_z_0_xxxyyz_xxzz[k] = -g_z_0_xxyyz_xxzz[k] * cd_x[k] + g_z_0_xxyyz_xxxzz[k];

                g_z_0_xxxyyz_xyyy[k] = -g_z_0_xxyyz_xyyy[k] * cd_x[k] + g_z_0_xxyyz_xxyyy[k];

                g_z_0_xxxyyz_xyyz[k] = -g_z_0_xxyyz_xyyz[k] * cd_x[k] + g_z_0_xxyyz_xxyyz[k];

                g_z_0_xxxyyz_xyzz[k] = -g_z_0_xxyyz_xyzz[k] * cd_x[k] + g_z_0_xxyyz_xxyzz[k];

                g_z_0_xxxyyz_xzzz[k] = -g_z_0_xxyyz_xzzz[k] * cd_x[k] + g_z_0_xxyyz_xxzzz[k];

                g_z_0_xxxyyz_yyyy[k] = -g_z_0_xxyyz_yyyy[k] * cd_x[k] + g_z_0_xxyyz_xyyyy[k];

                g_z_0_xxxyyz_yyyz[k] = -g_z_0_xxyyz_yyyz[k] * cd_x[k] + g_z_0_xxyyz_xyyyz[k];

                g_z_0_xxxyyz_yyzz[k] = -g_z_0_xxyyz_yyzz[k] * cd_x[k] + g_z_0_xxyyz_xyyzz[k];

                g_z_0_xxxyyz_yzzz[k] = -g_z_0_xxyyz_yzzz[k] * cd_x[k] + g_z_0_xxyyz_xyzzz[k];

                g_z_0_xxxyyz_zzzz[k] = -g_z_0_xxyyz_zzzz[k] * cd_x[k] + g_z_0_xxyyz_xzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 120);

            auto g_z_0_xxxyzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 121);

            auto g_z_0_xxxyzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 122);

            auto g_z_0_xxxyzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 123);

            auto g_z_0_xxxyzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 124);

            auto g_z_0_xxxyzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 125);

            auto g_z_0_xxxyzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 126);

            auto g_z_0_xxxyzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 127);

            auto g_z_0_xxxyzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 128);

            auto g_z_0_xxxyzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 129);

            auto g_z_0_xxxyzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 130);

            auto g_z_0_xxxyzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 131);

            auto g_z_0_xxxyzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 132);

            auto g_z_0_xxxyzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 133);

            auto g_z_0_xxxyzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 134);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyzz_xxxx, g_z_0_xxxyzz_xxxy, g_z_0_xxxyzz_xxxz, g_z_0_xxxyzz_xxyy, g_z_0_xxxyzz_xxyz, g_z_0_xxxyzz_xxzz, g_z_0_xxxyzz_xyyy, g_z_0_xxxyzz_xyyz, g_z_0_xxxyzz_xyzz, g_z_0_xxxyzz_xzzz, g_z_0_xxxyzz_yyyy, g_z_0_xxxyzz_yyyz, g_z_0_xxxyzz_yyzz, g_z_0_xxxyzz_yzzz, g_z_0_xxxyzz_zzzz, g_z_0_xxyzz_xxxx, g_z_0_xxyzz_xxxxx, g_z_0_xxyzz_xxxxy, g_z_0_xxyzz_xxxxz, g_z_0_xxyzz_xxxy, g_z_0_xxyzz_xxxyy, g_z_0_xxyzz_xxxyz, g_z_0_xxyzz_xxxz, g_z_0_xxyzz_xxxzz, g_z_0_xxyzz_xxyy, g_z_0_xxyzz_xxyyy, g_z_0_xxyzz_xxyyz, g_z_0_xxyzz_xxyz, g_z_0_xxyzz_xxyzz, g_z_0_xxyzz_xxzz, g_z_0_xxyzz_xxzzz, g_z_0_xxyzz_xyyy, g_z_0_xxyzz_xyyyy, g_z_0_xxyzz_xyyyz, g_z_0_xxyzz_xyyz, g_z_0_xxyzz_xyyzz, g_z_0_xxyzz_xyzz, g_z_0_xxyzz_xyzzz, g_z_0_xxyzz_xzzz, g_z_0_xxyzz_xzzzz, g_z_0_xxyzz_yyyy, g_z_0_xxyzz_yyyz, g_z_0_xxyzz_yyzz, g_z_0_xxyzz_yzzz, g_z_0_xxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_xxxx[k] = -g_z_0_xxyzz_xxxx[k] * cd_x[k] + g_z_0_xxyzz_xxxxx[k];

                g_z_0_xxxyzz_xxxy[k] = -g_z_0_xxyzz_xxxy[k] * cd_x[k] + g_z_0_xxyzz_xxxxy[k];

                g_z_0_xxxyzz_xxxz[k] = -g_z_0_xxyzz_xxxz[k] * cd_x[k] + g_z_0_xxyzz_xxxxz[k];

                g_z_0_xxxyzz_xxyy[k] = -g_z_0_xxyzz_xxyy[k] * cd_x[k] + g_z_0_xxyzz_xxxyy[k];

                g_z_0_xxxyzz_xxyz[k] = -g_z_0_xxyzz_xxyz[k] * cd_x[k] + g_z_0_xxyzz_xxxyz[k];

                g_z_0_xxxyzz_xxzz[k] = -g_z_0_xxyzz_xxzz[k] * cd_x[k] + g_z_0_xxyzz_xxxzz[k];

                g_z_0_xxxyzz_xyyy[k] = -g_z_0_xxyzz_xyyy[k] * cd_x[k] + g_z_0_xxyzz_xxyyy[k];

                g_z_0_xxxyzz_xyyz[k] = -g_z_0_xxyzz_xyyz[k] * cd_x[k] + g_z_0_xxyzz_xxyyz[k];

                g_z_0_xxxyzz_xyzz[k] = -g_z_0_xxyzz_xyzz[k] * cd_x[k] + g_z_0_xxyzz_xxyzz[k];

                g_z_0_xxxyzz_xzzz[k] = -g_z_0_xxyzz_xzzz[k] * cd_x[k] + g_z_0_xxyzz_xxzzz[k];

                g_z_0_xxxyzz_yyyy[k] = -g_z_0_xxyzz_yyyy[k] * cd_x[k] + g_z_0_xxyzz_xyyyy[k];

                g_z_0_xxxyzz_yyyz[k] = -g_z_0_xxyzz_yyyz[k] * cd_x[k] + g_z_0_xxyzz_xyyyz[k];

                g_z_0_xxxyzz_yyzz[k] = -g_z_0_xxyzz_yyzz[k] * cd_x[k] + g_z_0_xxyzz_xyyzz[k];

                g_z_0_xxxyzz_yzzz[k] = -g_z_0_xxyzz_yzzz[k] * cd_x[k] + g_z_0_xxyzz_xyzzz[k];

                g_z_0_xxxyzz_zzzz[k] = -g_z_0_xxyzz_zzzz[k] * cd_x[k] + g_z_0_xxyzz_xzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 135);

            auto g_z_0_xxxzzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 136);

            auto g_z_0_xxxzzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 137);

            auto g_z_0_xxxzzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 138);

            auto g_z_0_xxxzzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 139);

            auto g_z_0_xxxzzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 140);

            auto g_z_0_xxxzzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 141);

            auto g_z_0_xxxzzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 142);

            auto g_z_0_xxxzzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 143);

            auto g_z_0_xxxzzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 144);

            auto g_z_0_xxxzzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 145);

            auto g_z_0_xxxzzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 146);

            auto g_z_0_xxxzzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 147);

            auto g_z_0_xxxzzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 148);

            auto g_z_0_xxxzzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzzz_xxxx, g_z_0_xxxzzz_xxxy, g_z_0_xxxzzz_xxxz, g_z_0_xxxzzz_xxyy, g_z_0_xxxzzz_xxyz, g_z_0_xxxzzz_xxzz, g_z_0_xxxzzz_xyyy, g_z_0_xxxzzz_xyyz, g_z_0_xxxzzz_xyzz, g_z_0_xxxzzz_xzzz, g_z_0_xxxzzz_yyyy, g_z_0_xxxzzz_yyyz, g_z_0_xxxzzz_yyzz, g_z_0_xxxzzz_yzzz, g_z_0_xxxzzz_zzzz, g_z_0_xxzzz_xxxx, g_z_0_xxzzz_xxxxx, g_z_0_xxzzz_xxxxy, g_z_0_xxzzz_xxxxz, g_z_0_xxzzz_xxxy, g_z_0_xxzzz_xxxyy, g_z_0_xxzzz_xxxyz, g_z_0_xxzzz_xxxz, g_z_0_xxzzz_xxxzz, g_z_0_xxzzz_xxyy, g_z_0_xxzzz_xxyyy, g_z_0_xxzzz_xxyyz, g_z_0_xxzzz_xxyz, g_z_0_xxzzz_xxyzz, g_z_0_xxzzz_xxzz, g_z_0_xxzzz_xxzzz, g_z_0_xxzzz_xyyy, g_z_0_xxzzz_xyyyy, g_z_0_xxzzz_xyyyz, g_z_0_xxzzz_xyyz, g_z_0_xxzzz_xyyzz, g_z_0_xxzzz_xyzz, g_z_0_xxzzz_xyzzz, g_z_0_xxzzz_xzzz, g_z_0_xxzzz_xzzzz, g_z_0_xxzzz_yyyy, g_z_0_xxzzz_yyyz, g_z_0_xxzzz_yyzz, g_z_0_xxzzz_yzzz, g_z_0_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_xxxx[k] = -g_z_0_xxzzz_xxxx[k] * cd_x[k] + g_z_0_xxzzz_xxxxx[k];

                g_z_0_xxxzzz_xxxy[k] = -g_z_0_xxzzz_xxxy[k] * cd_x[k] + g_z_0_xxzzz_xxxxy[k];

                g_z_0_xxxzzz_xxxz[k] = -g_z_0_xxzzz_xxxz[k] * cd_x[k] + g_z_0_xxzzz_xxxxz[k];

                g_z_0_xxxzzz_xxyy[k] = -g_z_0_xxzzz_xxyy[k] * cd_x[k] + g_z_0_xxzzz_xxxyy[k];

                g_z_0_xxxzzz_xxyz[k] = -g_z_0_xxzzz_xxyz[k] * cd_x[k] + g_z_0_xxzzz_xxxyz[k];

                g_z_0_xxxzzz_xxzz[k] = -g_z_0_xxzzz_xxzz[k] * cd_x[k] + g_z_0_xxzzz_xxxzz[k];

                g_z_0_xxxzzz_xyyy[k] = -g_z_0_xxzzz_xyyy[k] * cd_x[k] + g_z_0_xxzzz_xxyyy[k];

                g_z_0_xxxzzz_xyyz[k] = -g_z_0_xxzzz_xyyz[k] * cd_x[k] + g_z_0_xxzzz_xxyyz[k];

                g_z_0_xxxzzz_xyzz[k] = -g_z_0_xxzzz_xyzz[k] * cd_x[k] + g_z_0_xxzzz_xxyzz[k];

                g_z_0_xxxzzz_xzzz[k] = -g_z_0_xxzzz_xzzz[k] * cd_x[k] + g_z_0_xxzzz_xxzzz[k];

                g_z_0_xxxzzz_yyyy[k] = -g_z_0_xxzzz_yyyy[k] * cd_x[k] + g_z_0_xxzzz_xyyyy[k];

                g_z_0_xxxzzz_yyyz[k] = -g_z_0_xxzzz_yyyz[k] * cd_x[k] + g_z_0_xxzzz_xyyyz[k];

                g_z_0_xxxzzz_yyzz[k] = -g_z_0_xxzzz_yyzz[k] * cd_x[k] + g_z_0_xxzzz_xyyzz[k];

                g_z_0_xxxzzz_yzzz[k] = -g_z_0_xxzzz_yzzz[k] * cd_x[k] + g_z_0_xxzzz_xyzzz[k];

                g_z_0_xxxzzz_zzzz[k] = -g_z_0_xxzzz_zzzz[k] * cd_x[k] + g_z_0_xxzzz_xzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 150);

            auto g_z_0_xxyyyy_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 151);

            auto g_z_0_xxyyyy_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 152);

            auto g_z_0_xxyyyy_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 153);

            auto g_z_0_xxyyyy_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 154);

            auto g_z_0_xxyyyy_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 155);

            auto g_z_0_xxyyyy_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 156);

            auto g_z_0_xxyyyy_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 157);

            auto g_z_0_xxyyyy_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 158);

            auto g_z_0_xxyyyy_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 159);

            auto g_z_0_xxyyyy_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 160);

            auto g_z_0_xxyyyy_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 161);

            auto g_z_0_xxyyyy_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 162);

            auto g_z_0_xxyyyy_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 163);

            auto g_z_0_xxyyyy_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 164);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyy_xxxx, g_z_0_xxyyyy_xxxy, g_z_0_xxyyyy_xxxz, g_z_0_xxyyyy_xxyy, g_z_0_xxyyyy_xxyz, g_z_0_xxyyyy_xxzz, g_z_0_xxyyyy_xyyy, g_z_0_xxyyyy_xyyz, g_z_0_xxyyyy_xyzz, g_z_0_xxyyyy_xzzz, g_z_0_xxyyyy_yyyy, g_z_0_xxyyyy_yyyz, g_z_0_xxyyyy_yyzz, g_z_0_xxyyyy_yzzz, g_z_0_xxyyyy_zzzz, g_z_0_xyyyy_xxxx, g_z_0_xyyyy_xxxxx, g_z_0_xyyyy_xxxxy, g_z_0_xyyyy_xxxxz, g_z_0_xyyyy_xxxy, g_z_0_xyyyy_xxxyy, g_z_0_xyyyy_xxxyz, g_z_0_xyyyy_xxxz, g_z_0_xyyyy_xxxzz, g_z_0_xyyyy_xxyy, g_z_0_xyyyy_xxyyy, g_z_0_xyyyy_xxyyz, g_z_0_xyyyy_xxyz, g_z_0_xyyyy_xxyzz, g_z_0_xyyyy_xxzz, g_z_0_xyyyy_xxzzz, g_z_0_xyyyy_xyyy, g_z_0_xyyyy_xyyyy, g_z_0_xyyyy_xyyyz, g_z_0_xyyyy_xyyz, g_z_0_xyyyy_xyyzz, g_z_0_xyyyy_xyzz, g_z_0_xyyyy_xyzzz, g_z_0_xyyyy_xzzz, g_z_0_xyyyy_xzzzz, g_z_0_xyyyy_yyyy, g_z_0_xyyyy_yyyz, g_z_0_xyyyy_yyzz, g_z_0_xyyyy_yzzz, g_z_0_xyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_xxxx[k] = -g_z_0_xyyyy_xxxx[k] * cd_x[k] + g_z_0_xyyyy_xxxxx[k];

                g_z_0_xxyyyy_xxxy[k] = -g_z_0_xyyyy_xxxy[k] * cd_x[k] + g_z_0_xyyyy_xxxxy[k];

                g_z_0_xxyyyy_xxxz[k] = -g_z_0_xyyyy_xxxz[k] * cd_x[k] + g_z_0_xyyyy_xxxxz[k];

                g_z_0_xxyyyy_xxyy[k] = -g_z_0_xyyyy_xxyy[k] * cd_x[k] + g_z_0_xyyyy_xxxyy[k];

                g_z_0_xxyyyy_xxyz[k] = -g_z_0_xyyyy_xxyz[k] * cd_x[k] + g_z_0_xyyyy_xxxyz[k];

                g_z_0_xxyyyy_xxzz[k] = -g_z_0_xyyyy_xxzz[k] * cd_x[k] + g_z_0_xyyyy_xxxzz[k];

                g_z_0_xxyyyy_xyyy[k] = -g_z_0_xyyyy_xyyy[k] * cd_x[k] + g_z_0_xyyyy_xxyyy[k];

                g_z_0_xxyyyy_xyyz[k] = -g_z_0_xyyyy_xyyz[k] * cd_x[k] + g_z_0_xyyyy_xxyyz[k];

                g_z_0_xxyyyy_xyzz[k] = -g_z_0_xyyyy_xyzz[k] * cd_x[k] + g_z_0_xyyyy_xxyzz[k];

                g_z_0_xxyyyy_xzzz[k] = -g_z_0_xyyyy_xzzz[k] * cd_x[k] + g_z_0_xyyyy_xxzzz[k];

                g_z_0_xxyyyy_yyyy[k] = -g_z_0_xyyyy_yyyy[k] * cd_x[k] + g_z_0_xyyyy_xyyyy[k];

                g_z_0_xxyyyy_yyyz[k] = -g_z_0_xyyyy_yyyz[k] * cd_x[k] + g_z_0_xyyyy_xyyyz[k];

                g_z_0_xxyyyy_yyzz[k] = -g_z_0_xyyyy_yyzz[k] * cd_x[k] + g_z_0_xyyyy_xyyzz[k];

                g_z_0_xxyyyy_yzzz[k] = -g_z_0_xyyyy_yzzz[k] * cd_x[k] + g_z_0_xyyyy_xyzzz[k];

                g_z_0_xxyyyy_zzzz[k] = -g_z_0_xyyyy_zzzz[k] * cd_x[k] + g_z_0_xyyyy_xzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 165);

            auto g_z_0_xxyyyz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 166);

            auto g_z_0_xxyyyz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 167);

            auto g_z_0_xxyyyz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 168);

            auto g_z_0_xxyyyz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 169);

            auto g_z_0_xxyyyz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 170);

            auto g_z_0_xxyyyz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 171);

            auto g_z_0_xxyyyz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 172);

            auto g_z_0_xxyyyz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 173);

            auto g_z_0_xxyyyz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 174);

            auto g_z_0_xxyyyz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 175);

            auto g_z_0_xxyyyz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 176);

            auto g_z_0_xxyyyz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 177);

            auto g_z_0_xxyyyz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 178);

            auto g_z_0_xxyyyz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyz_xxxx, g_z_0_xxyyyz_xxxy, g_z_0_xxyyyz_xxxz, g_z_0_xxyyyz_xxyy, g_z_0_xxyyyz_xxyz, g_z_0_xxyyyz_xxzz, g_z_0_xxyyyz_xyyy, g_z_0_xxyyyz_xyyz, g_z_0_xxyyyz_xyzz, g_z_0_xxyyyz_xzzz, g_z_0_xxyyyz_yyyy, g_z_0_xxyyyz_yyyz, g_z_0_xxyyyz_yyzz, g_z_0_xxyyyz_yzzz, g_z_0_xxyyyz_zzzz, g_z_0_xyyyz_xxxx, g_z_0_xyyyz_xxxxx, g_z_0_xyyyz_xxxxy, g_z_0_xyyyz_xxxxz, g_z_0_xyyyz_xxxy, g_z_0_xyyyz_xxxyy, g_z_0_xyyyz_xxxyz, g_z_0_xyyyz_xxxz, g_z_0_xyyyz_xxxzz, g_z_0_xyyyz_xxyy, g_z_0_xyyyz_xxyyy, g_z_0_xyyyz_xxyyz, g_z_0_xyyyz_xxyz, g_z_0_xyyyz_xxyzz, g_z_0_xyyyz_xxzz, g_z_0_xyyyz_xxzzz, g_z_0_xyyyz_xyyy, g_z_0_xyyyz_xyyyy, g_z_0_xyyyz_xyyyz, g_z_0_xyyyz_xyyz, g_z_0_xyyyz_xyyzz, g_z_0_xyyyz_xyzz, g_z_0_xyyyz_xyzzz, g_z_0_xyyyz_xzzz, g_z_0_xyyyz_xzzzz, g_z_0_xyyyz_yyyy, g_z_0_xyyyz_yyyz, g_z_0_xyyyz_yyzz, g_z_0_xyyyz_yzzz, g_z_0_xyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_xxxx[k] = -g_z_0_xyyyz_xxxx[k] * cd_x[k] + g_z_0_xyyyz_xxxxx[k];

                g_z_0_xxyyyz_xxxy[k] = -g_z_0_xyyyz_xxxy[k] * cd_x[k] + g_z_0_xyyyz_xxxxy[k];

                g_z_0_xxyyyz_xxxz[k] = -g_z_0_xyyyz_xxxz[k] * cd_x[k] + g_z_0_xyyyz_xxxxz[k];

                g_z_0_xxyyyz_xxyy[k] = -g_z_0_xyyyz_xxyy[k] * cd_x[k] + g_z_0_xyyyz_xxxyy[k];

                g_z_0_xxyyyz_xxyz[k] = -g_z_0_xyyyz_xxyz[k] * cd_x[k] + g_z_0_xyyyz_xxxyz[k];

                g_z_0_xxyyyz_xxzz[k] = -g_z_0_xyyyz_xxzz[k] * cd_x[k] + g_z_0_xyyyz_xxxzz[k];

                g_z_0_xxyyyz_xyyy[k] = -g_z_0_xyyyz_xyyy[k] * cd_x[k] + g_z_0_xyyyz_xxyyy[k];

                g_z_0_xxyyyz_xyyz[k] = -g_z_0_xyyyz_xyyz[k] * cd_x[k] + g_z_0_xyyyz_xxyyz[k];

                g_z_0_xxyyyz_xyzz[k] = -g_z_0_xyyyz_xyzz[k] * cd_x[k] + g_z_0_xyyyz_xxyzz[k];

                g_z_0_xxyyyz_xzzz[k] = -g_z_0_xyyyz_xzzz[k] * cd_x[k] + g_z_0_xyyyz_xxzzz[k];

                g_z_0_xxyyyz_yyyy[k] = -g_z_0_xyyyz_yyyy[k] * cd_x[k] + g_z_0_xyyyz_xyyyy[k];

                g_z_0_xxyyyz_yyyz[k] = -g_z_0_xyyyz_yyyz[k] * cd_x[k] + g_z_0_xyyyz_xyyyz[k];

                g_z_0_xxyyyz_yyzz[k] = -g_z_0_xyyyz_yyzz[k] * cd_x[k] + g_z_0_xyyyz_xyyzz[k];

                g_z_0_xxyyyz_yzzz[k] = -g_z_0_xyyyz_yzzz[k] * cd_x[k] + g_z_0_xyyyz_xyzzz[k];

                g_z_0_xxyyyz_zzzz[k] = -g_z_0_xyyyz_zzzz[k] * cd_x[k] + g_z_0_xyyyz_xzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 180);

            auto g_z_0_xxyyzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 181);

            auto g_z_0_xxyyzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 182);

            auto g_z_0_xxyyzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 183);

            auto g_z_0_xxyyzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 184);

            auto g_z_0_xxyyzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 185);

            auto g_z_0_xxyyzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 186);

            auto g_z_0_xxyyzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 187);

            auto g_z_0_xxyyzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 188);

            auto g_z_0_xxyyzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 189);

            auto g_z_0_xxyyzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 190);

            auto g_z_0_xxyyzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 191);

            auto g_z_0_xxyyzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 192);

            auto g_z_0_xxyyzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 193);

            auto g_z_0_xxyyzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 194);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyzz_xxxx, g_z_0_xxyyzz_xxxy, g_z_0_xxyyzz_xxxz, g_z_0_xxyyzz_xxyy, g_z_0_xxyyzz_xxyz, g_z_0_xxyyzz_xxzz, g_z_0_xxyyzz_xyyy, g_z_0_xxyyzz_xyyz, g_z_0_xxyyzz_xyzz, g_z_0_xxyyzz_xzzz, g_z_0_xxyyzz_yyyy, g_z_0_xxyyzz_yyyz, g_z_0_xxyyzz_yyzz, g_z_0_xxyyzz_yzzz, g_z_0_xxyyzz_zzzz, g_z_0_xyyzz_xxxx, g_z_0_xyyzz_xxxxx, g_z_0_xyyzz_xxxxy, g_z_0_xyyzz_xxxxz, g_z_0_xyyzz_xxxy, g_z_0_xyyzz_xxxyy, g_z_0_xyyzz_xxxyz, g_z_0_xyyzz_xxxz, g_z_0_xyyzz_xxxzz, g_z_0_xyyzz_xxyy, g_z_0_xyyzz_xxyyy, g_z_0_xyyzz_xxyyz, g_z_0_xyyzz_xxyz, g_z_0_xyyzz_xxyzz, g_z_0_xyyzz_xxzz, g_z_0_xyyzz_xxzzz, g_z_0_xyyzz_xyyy, g_z_0_xyyzz_xyyyy, g_z_0_xyyzz_xyyyz, g_z_0_xyyzz_xyyz, g_z_0_xyyzz_xyyzz, g_z_0_xyyzz_xyzz, g_z_0_xyyzz_xyzzz, g_z_0_xyyzz_xzzz, g_z_0_xyyzz_xzzzz, g_z_0_xyyzz_yyyy, g_z_0_xyyzz_yyyz, g_z_0_xyyzz_yyzz, g_z_0_xyyzz_yzzz, g_z_0_xyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_xxxx[k] = -g_z_0_xyyzz_xxxx[k] * cd_x[k] + g_z_0_xyyzz_xxxxx[k];

                g_z_0_xxyyzz_xxxy[k] = -g_z_0_xyyzz_xxxy[k] * cd_x[k] + g_z_0_xyyzz_xxxxy[k];

                g_z_0_xxyyzz_xxxz[k] = -g_z_0_xyyzz_xxxz[k] * cd_x[k] + g_z_0_xyyzz_xxxxz[k];

                g_z_0_xxyyzz_xxyy[k] = -g_z_0_xyyzz_xxyy[k] * cd_x[k] + g_z_0_xyyzz_xxxyy[k];

                g_z_0_xxyyzz_xxyz[k] = -g_z_0_xyyzz_xxyz[k] * cd_x[k] + g_z_0_xyyzz_xxxyz[k];

                g_z_0_xxyyzz_xxzz[k] = -g_z_0_xyyzz_xxzz[k] * cd_x[k] + g_z_0_xyyzz_xxxzz[k];

                g_z_0_xxyyzz_xyyy[k] = -g_z_0_xyyzz_xyyy[k] * cd_x[k] + g_z_0_xyyzz_xxyyy[k];

                g_z_0_xxyyzz_xyyz[k] = -g_z_0_xyyzz_xyyz[k] * cd_x[k] + g_z_0_xyyzz_xxyyz[k];

                g_z_0_xxyyzz_xyzz[k] = -g_z_0_xyyzz_xyzz[k] * cd_x[k] + g_z_0_xyyzz_xxyzz[k];

                g_z_0_xxyyzz_xzzz[k] = -g_z_0_xyyzz_xzzz[k] * cd_x[k] + g_z_0_xyyzz_xxzzz[k];

                g_z_0_xxyyzz_yyyy[k] = -g_z_0_xyyzz_yyyy[k] * cd_x[k] + g_z_0_xyyzz_xyyyy[k];

                g_z_0_xxyyzz_yyyz[k] = -g_z_0_xyyzz_yyyz[k] * cd_x[k] + g_z_0_xyyzz_xyyyz[k];

                g_z_0_xxyyzz_yyzz[k] = -g_z_0_xyyzz_yyzz[k] * cd_x[k] + g_z_0_xyyzz_xyyzz[k];

                g_z_0_xxyyzz_yzzz[k] = -g_z_0_xyyzz_yzzz[k] * cd_x[k] + g_z_0_xyyzz_xyzzz[k];

                g_z_0_xxyyzz_zzzz[k] = -g_z_0_xyyzz_zzzz[k] * cd_x[k] + g_z_0_xyyzz_xzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 195);

            auto g_z_0_xxyzzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 196);

            auto g_z_0_xxyzzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 197);

            auto g_z_0_xxyzzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 198);

            auto g_z_0_xxyzzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 199);

            auto g_z_0_xxyzzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 200);

            auto g_z_0_xxyzzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 201);

            auto g_z_0_xxyzzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 202);

            auto g_z_0_xxyzzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 203);

            auto g_z_0_xxyzzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 204);

            auto g_z_0_xxyzzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 205);

            auto g_z_0_xxyzzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 206);

            auto g_z_0_xxyzzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 207);

            auto g_z_0_xxyzzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 208);

            auto g_z_0_xxyzzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzzz_xxxx, g_z_0_xxyzzz_xxxy, g_z_0_xxyzzz_xxxz, g_z_0_xxyzzz_xxyy, g_z_0_xxyzzz_xxyz, g_z_0_xxyzzz_xxzz, g_z_0_xxyzzz_xyyy, g_z_0_xxyzzz_xyyz, g_z_0_xxyzzz_xyzz, g_z_0_xxyzzz_xzzz, g_z_0_xxyzzz_yyyy, g_z_0_xxyzzz_yyyz, g_z_0_xxyzzz_yyzz, g_z_0_xxyzzz_yzzz, g_z_0_xxyzzz_zzzz, g_z_0_xyzzz_xxxx, g_z_0_xyzzz_xxxxx, g_z_0_xyzzz_xxxxy, g_z_0_xyzzz_xxxxz, g_z_0_xyzzz_xxxy, g_z_0_xyzzz_xxxyy, g_z_0_xyzzz_xxxyz, g_z_0_xyzzz_xxxz, g_z_0_xyzzz_xxxzz, g_z_0_xyzzz_xxyy, g_z_0_xyzzz_xxyyy, g_z_0_xyzzz_xxyyz, g_z_0_xyzzz_xxyz, g_z_0_xyzzz_xxyzz, g_z_0_xyzzz_xxzz, g_z_0_xyzzz_xxzzz, g_z_0_xyzzz_xyyy, g_z_0_xyzzz_xyyyy, g_z_0_xyzzz_xyyyz, g_z_0_xyzzz_xyyz, g_z_0_xyzzz_xyyzz, g_z_0_xyzzz_xyzz, g_z_0_xyzzz_xyzzz, g_z_0_xyzzz_xzzz, g_z_0_xyzzz_xzzzz, g_z_0_xyzzz_yyyy, g_z_0_xyzzz_yyyz, g_z_0_xyzzz_yyzz, g_z_0_xyzzz_yzzz, g_z_0_xyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_xxxx[k] = -g_z_0_xyzzz_xxxx[k] * cd_x[k] + g_z_0_xyzzz_xxxxx[k];

                g_z_0_xxyzzz_xxxy[k] = -g_z_0_xyzzz_xxxy[k] * cd_x[k] + g_z_0_xyzzz_xxxxy[k];

                g_z_0_xxyzzz_xxxz[k] = -g_z_0_xyzzz_xxxz[k] * cd_x[k] + g_z_0_xyzzz_xxxxz[k];

                g_z_0_xxyzzz_xxyy[k] = -g_z_0_xyzzz_xxyy[k] * cd_x[k] + g_z_0_xyzzz_xxxyy[k];

                g_z_0_xxyzzz_xxyz[k] = -g_z_0_xyzzz_xxyz[k] * cd_x[k] + g_z_0_xyzzz_xxxyz[k];

                g_z_0_xxyzzz_xxzz[k] = -g_z_0_xyzzz_xxzz[k] * cd_x[k] + g_z_0_xyzzz_xxxzz[k];

                g_z_0_xxyzzz_xyyy[k] = -g_z_0_xyzzz_xyyy[k] * cd_x[k] + g_z_0_xyzzz_xxyyy[k];

                g_z_0_xxyzzz_xyyz[k] = -g_z_0_xyzzz_xyyz[k] * cd_x[k] + g_z_0_xyzzz_xxyyz[k];

                g_z_0_xxyzzz_xyzz[k] = -g_z_0_xyzzz_xyzz[k] * cd_x[k] + g_z_0_xyzzz_xxyzz[k];

                g_z_0_xxyzzz_xzzz[k] = -g_z_0_xyzzz_xzzz[k] * cd_x[k] + g_z_0_xyzzz_xxzzz[k];

                g_z_0_xxyzzz_yyyy[k] = -g_z_0_xyzzz_yyyy[k] * cd_x[k] + g_z_0_xyzzz_xyyyy[k];

                g_z_0_xxyzzz_yyyz[k] = -g_z_0_xyzzz_yyyz[k] * cd_x[k] + g_z_0_xyzzz_xyyyz[k];

                g_z_0_xxyzzz_yyzz[k] = -g_z_0_xyzzz_yyzz[k] * cd_x[k] + g_z_0_xyzzz_xyyzz[k];

                g_z_0_xxyzzz_yzzz[k] = -g_z_0_xyzzz_yzzz[k] * cd_x[k] + g_z_0_xyzzz_xyzzz[k];

                g_z_0_xxyzzz_zzzz[k] = -g_z_0_xyzzz_zzzz[k] * cd_x[k] + g_z_0_xyzzz_xzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 210);

            auto g_z_0_xxzzzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 211);

            auto g_z_0_xxzzzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 212);

            auto g_z_0_xxzzzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 213);

            auto g_z_0_xxzzzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 214);

            auto g_z_0_xxzzzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 215);

            auto g_z_0_xxzzzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 216);

            auto g_z_0_xxzzzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 217);

            auto g_z_0_xxzzzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 218);

            auto g_z_0_xxzzzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 219);

            auto g_z_0_xxzzzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 220);

            auto g_z_0_xxzzzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 221);

            auto g_z_0_xxzzzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 222);

            auto g_z_0_xxzzzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 223);

            auto g_z_0_xxzzzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 224);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzzz_xxxx, g_z_0_xxzzzz_xxxy, g_z_0_xxzzzz_xxxz, g_z_0_xxzzzz_xxyy, g_z_0_xxzzzz_xxyz, g_z_0_xxzzzz_xxzz, g_z_0_xxzzzz_xyyy, g_z_0_xxzzzz_xyyz, g_z_0_xxzzzz_xyzz, g_z_0_xxzzzz_xzzz, g_z_0_xxzzzz_yyyy, g_z_0_xxzzzz_yyyz, g_z_0_xxzzzz_yyzz, g_z_0_xxzzzz_yzzz, g_z_0_xxzzzz_zzzz, g_z_0_xzzzz_xxxx, g_z_0_xzzzz_xxxxx, g_z_0_xzzzz_xxxxy, g_z_0_xzzzz_xxxxz, g_z_0_xzzzz_xxxy, g_z_0_xzzzz_xxxyy, g_z_0_xzzzz_xxxyz, g_z_0_xzzzz_xxxz, g_z_0_xzzzz_xxxzz, g_z_0_xzzzz_xxyy, g_z_0_xzzzz_xxyyy, g_z_0_xzzzz_xxyyz, g_z_0_xzzzz_xxyz, g_z_0_xzzzz_xxyzz, g_z_0_xzzzz_xxzz, g_z_0_xzzzz_xxzzz, g_z_0_xzzzz_xyyy, g_z_0_xzzzz_xyyyy, g_z_0_xzzzz_xyyyz, g_z_0_xzzzz_xyyz, g_z_0_xzzzz_xyyzz, g_z_0_xzzzz_xyzz, g_z_0_xzzzz_xyzzz, g_z_0_xzzzz_xzzz, g_z_0_xzzzz_xzzzz, g_z_0_xzzzz_yyyy, g_z_0_xzzzz_yyyz, g_z_0_xzzzz_yyzz, g_z_0_xzzzz_yzzz, g_z_0_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_xxxx[k] = -g_z_0_xzzzz_xxxx[k] * cd_x[k] + g_z_0_xzzzz_xxxxx[k];

                g_z_0_xxzzzz_xxxy[k] = -g_z_0_xzzzz_xxxy[k] * cd_x[k] + g_z_0_xzzzz_xxxxy[k];

                g_z_0_xxzzzz_xxxz[k] = -g_z_0_xzzzz_xxxz[k] * cd_x[k] + g_z_0_xzzzz_xxxxz[k];

                g_z_0_xxzzzz_xxyy[k] = -g_z_0_xzzzz_xxyy[k] * cd_x[k] + g_z_0_xzzzz_xxxyy[k];

                g_z_0_xxzzzz_xxyz[k] = -g_z_0_xzzzz_xxyz[k] * cd_x[k] + g_z_0_xzzzz_xxxyz[k];

                g_z_0_xxzzzz_xxzz[k] = -g_z_0_xzzzz_xxzz[k] * cd_x[k] + g_z_0_xzzzz_xxxzz[k];

                g_z_0_xxzzzz_xyyy[k] = -g_z_0_xzzzz_xyyy[k] * cd_x[k] + g_z_0_xzzzz_xxyyy[k];

                g_z_0_xxzzzz_xyyz[k] = -g_z_0_xzzzz_xyyz[k] * cd_x[k] + g_z_0_xzzzz_xxyyz[k];

                g_z_0_xxzzzz_xyzz[k] = -g_z_0_xzzzz_xyzz[k] * cd_x[k] + g_z_0_xzzzz_xxyzz[k];

                g_z_0_xxzzzz_xzzz[k] = -g_z_0_xzzzz_xzzz[k] * cd_x[k] + g_z_0_xzzzz_xxzzz[k];

                g_z_0_xxzzzz_yyyy[k] = -g_z_0_xzzzz_yyyy[k] * cd_x[k] + g_z_0_xzzzz_xyyyy[k];

                g_z_0_xxzzzz_yyyz[k] = -g_z_0_xzzzz_yyyz[k] * cd_x[k] + g_z_0_xzzzz_xyyyz[k];

                g_z_0_xxzzzz_yyzz[k] = -g_z_0_xzzzz_yyzz[k] * cd_x[k] + g_z_0_xzzzz_xyyzz[k];

                g_z_0_xxzzzz_yzzz[k] = -g_z_0_xzzzz_yzzz[k] * cd_x[k] + g_z_0_xzzzz_xyzzz[k];

                g_z_0_xxzzzz_zzzz[k] = -g_z_0_xzzzz_zzzz[k] * cd_x[k] + g_z_0_xzzzz_xzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 225);

            auto g_z_0_xyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 226);

            auto g_z_0_xyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 227);

            auto g_z_0_xyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 228);

            auto g_z_0_xyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 229);

            auto g_z_0_xyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 230);

            auto g_z_0_xyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 231);

            auto g_z_0_xyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 232);

            auto g_z_0_xyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 233);

            auto g_z_0_xyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 234);

            auto g_z_0_xyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 235);

            auto g_z_0_xyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 236);

            auto g_z_0_xyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 237);

            auto g_z_0_xyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 238);

            auto g_z_0_xyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 239);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyy_xxxx, g_z_0_xyyyyy_xxxy, g_z_0_xyyyyy_xxxz, g_z_0_xyyyyy_xxyy, g_z_0_xyyyyy_xxyz, g_z_0_xyyyyy_xxzz, g_z_0_xyyyyy_xyyy, g_z_0_xyyyyy_xyyz, g_z_0_xyyyyy_xyzz, g_z_0_xyyyyy_xzzz, g_z_0_xyyyyy_yyyy, g_z_0_xyyyyy_yyyz, g_z_0_xyyyyy_yyzz, g_z_0_xyyyyy_yzzz, g_z_0_xyyyyy_zzzz, g_z_0_yyyyy_xxxx, g_z_0_yyyyy_xxxxx, g_z_0_yyyyy_xxxxy, g_z_0_yyyyy_xxxxz, g_z_0_yyyyy_xxxy, g_z_0_yyyyy_xxxyy, g_z_0_yyyyy_xxxyz, g_z_0_yyyyy_xxxz, g_z_0_yyyyy_xxxzz, g_z_0_yyyyy_xxyy, g_z_0_yyyyy_xxyyy, g_z_0_yyyyy_xxyyz, g_z_0_yyyyy_xxyz, g_z_0_yyyyy_xxyzz, g_z_0_yyyyy_xxzz, g_z_0_yyyyy_xxzzz, g_z_0_yyyyy_xyyy, g_z_0_yyyyy_xyyyy, g_z_0_yyyyy_xyyyz, g_z_0_yyyyy_xyyz, g_z_0_yyyyy_xyyzz, g_z_0_yyyyy_xyzz, g_z_0_yyyyy_xyzzz, g_z_0_yyyyy_xzzz, g_z_0_yyyyy_xzzzz, g_z_0_yyyyy_yyyy, g_z_0_yyyyy_yyyz, g_z_0_yyyyy_yyzz, g_z_0_yyyyy_yzzz, g_z_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_xxxx[k] = -g_z_0_yyyyy_xxxx[k] * cd_x[k] + g_z_0_yyyyy_xxxxx[k];

                g_z_0_xyyyyy_xxxy[k] = -g_z_0_yyyyy_xxxy[k] * cd_x[k] + g_z_0_yyyyy_xxxxy[k];

                g_z_0_xyyyyy_xxxz[k] = -g_z_0_yyyyy_xxxz[k] * cd_x[k] + g_z_0_yyyyy_xxxxz[k];

                g_z_0_xyyyyy_xxyy[k] = -g_z_0_yyyyy_xxyy[k] * cd_x[k] + g_z_0_yyyyy_xxxyy[k];

                g_z_0_xyyyyy_xxyz[k] = -g_z_0_yyyyy_xxyz[k] * cd_x[k] + g_z_0_yyyyy_xxxyz[k];

                g_z_0_xyyyyy_xxzz[k] = -g_z_0_yyyyy_xxzz[k] * cd_x[k] + g_z_0_yyyyy_xxxzz[k];

                g_z_0_xyyyyy_xyyy[k] = -g_z_0_yyyyy_xyyy[k] * cd_x[k] + g_z_0_yyyyy_xxyyy[k];

                g_z_0_xyyyyy_xyyz[k] = -g_z_0_yyyyy_xyyz[k] * cd_x[k] + g_z_0_yyyyy_xxyyz[k];

                g_z_0_xyyyyy_xyzz[k] = -g_z_0_yyyyy_xyzz[k] * cd_x[k] + g_z_0_yyyyy_xxyzz[k];

                g_z_0_xyyyyy_xzzz[k] = -g_z_0_yyyyy_xzzz[k] * cd_x[k] + g_z_0_yyyyy_xxzzz[k];

                g_z_0_xyyyyy_yyyy[k] = -g_z_0_yyyyy_yyyy[k] * cd_x[k] + g_z_0_yyyyy_xyyyy[k];

                g_z_0_xyyyyy_yyyz[k] = -g_z_0_yyyyy_yyyz[k] * cd_x[k] + g_z_0_yyyyy_xyyyz[k];

                g_z_0_xyyyyy_yyzz[k] = -g_z_0_yyyyy_yyzz[k] * cd_x[k] + g_z_0_yyyyy_xyyzz[k];

                g_z_0_xyyyyy_yzzz[k] = -g_z_0_yyyyy_yzzz[k] * cd_x[k] + g_z_0_yyyyy_xyzzz[k];

                g_z_0_xyyyyy_zzzz[k] = -g_z_0_yyyyy_zzzz[k] * cd_x[k] + g_z_0_yyyyy_xzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 240);

            auto g_z_0_xyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 241);

            auto g_z_0_xyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 242);

            auto g_z_0_xyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 243);

            auto g_z_0_xyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 244);

            auto g_z_0_xyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 245);

            auto g_z_0_xyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 246);

            auto g_z_0_xyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 247);

            auto g_z_0_xyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 248);

            auto g_z_0_xyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 249);

            auto g_z_0_xyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 250);

            auto g_z_0_xyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 251);

            auto g_z_0_xyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 252);

            auto g_z_0_xyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 253);

            auto g_z_0_xyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 254);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyz_xxxx, g_z_0_xyyyyz_xxxy, g_z_0_xyyyyz_xxxz, g_z_0_xyyyyz_xxyy, g_z_0_xyyyyz_xxyz, g_z_0_xyyyyz_xxzz, g_z_0_xyyyyz_xyyy, g_z_0_xyyyyz_xyyz, g_z_0_xyyyyz_xyzz, g_z_0_xyyyyz_xzzz, g_z_0_xyyyyz_yyyy, g_z_0_xyyyyz_yyyz, g_z_0_xyyyyz_yyzz, g_z_0_xyyyyz_yzzz, g_z_0_xyyyyz_zzzz, g_z_0_yyyyz_xxxx, g_z_0_yyyyz_xxxxx, g_z_0_yyyyz_xxxxy, g_z_0_yyyyz_xxxxz, g_z_0_yyyyz_xxxy, g_z_0_yyyyz_xxxyy, g_z_0_yyyyz_xxxyz, g_z_0_yyyyz_xxxz, g_z_0_yyyyz_xxxzz, g_z_0_yyyyz_xxyy, g_z_0_yyyyz_xxyyy, g_z_0_yyyyz_xxyyz, g_z_0_yyyyz_xxyz, g_z_0_yyyyz_xxyzz, g_z_0_yyyyz_xxzz, g_z_0_yyyyz_xxzzz, g_z_0_yyyyz_xyyy, g_z_0_yyyyz_xyyyy, g_z_0_yyyyz_xyyyz, g_z_0_yyyyz_xyyz, g_z_0_yyyyz_xyyzz, g_z_0_yyyyz_xyzz, g_z_0_yyyyz_xyzzz, g_z_0_yyyyz_xzzz, g_z_0_yyyyz_xzzzz, g_z_0_yyyyz_yyyy, g_z_0_yyyyz_yyyz, g_z_0_yyyyz_yyzz, g_z_0_yyyyz_yzzz, g_z_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_xxxx[k] = -g_z_0_yyyyz_xxxx[k] * cd_x[k] + g_z_0_yyyyz_xxxxx[k];

                g_z_0_xyyyyz_xxxy[k] = -g_z_0_yyyyz_xxxy[k] * cd_x[k] + g_z_0_yyyyz_xxxxy[k];

                g_z_0_xyyyyz_xxxz[k] = -g_z_0_yyyyz_xxxz[k] * cd_x[k] + g_z_0_yyyyz_xxxxz[k];

                g_z_0_xyyyyz_xxyy[k] = -g_z_0_yyyyz_xxyy[k] * cd_x[k] + g_z_0_yyyyz_xxxyy[k];

                g_z_0_xyyyyz_xxyz[k] = -g_z_0_yyyyz_xxyz[k] * cd_x[k] + g_z_0_yyyyz_xxxyz[k];

                g_z_0_xyyyyz_xxzz[k] = -g_z_0_yyyyz_xxzz[k] * cd_x[k] + g_z_0_yyyyz_xxxzz[k];

                g_z_0_xyyyyz_xyyy[k] = -g_z_0_yyyyz_xyyy[k] * cd_x[k] + g_z_0_yyyyz_xxyyy[k];

                g_z_0_xyyyyz_xyyz[k] = -g_z_0_yyyyz_xyyz[k] * cd_x[k] + g_z_0_yyyyz_xxyyz[k];

                g_z_0_xyyyyz_xyzz[k] = -g_z_0_yyyyz_xyzz[k] * cd_x[k] + g_z_0_yyyyz_xxyzz[k];

                g_z_0_xyyyyz_xzzz[k] = -g_z_0_yyyyz_xzzz[k] * cd_x[k] + g_z_0_yyyyz_xxzzz[k];

                g_z_0_xyyyyz_yyyy[k] = -g_z_0_yyyyz_yyyy[k] * cd_x[k] + g_z_0_yyyyz_xyyyy[k];

                g_z_0_xyyyyz_yyyz[k] = -g_z_0_yyyyz_yyyz[k] * cd_x[k] + g_z_0_yyyyz_xyyyz[k];

                g_z_0_xyyyyz_yyzz[k] = -g_z_0_yyyyz_yyzz[k] * cd_x[k] + g_z_0_yyyyz_xyyzz[k];

                g_z_0_xyyyyz_yzzz[k] = -g_z_0_yyyyz_yzzz[k] * cd_x[k] + g_z_0_yyyyz_xyzzz[k];

                g_z_0_xyyyyz_zzzz[k] = -g_z_0_yyyyz_zzzz[k] * cd_x[k] + g_z_0_yyyyz_xzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 255);

            auto g_z_0_xyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 256);

            auto g_z_0_xyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 257);

            auto g_z_0_xyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 258);

            auto g_z_0_xyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 259);

            auto g_z_0_xyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 260);

            auto g_z_0_xyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 261);

            auto g_z_0_xyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 262);

            auto g_z_0_xyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 263);

            auto g_z_0_xyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 264);

            auto g_z_0_xyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 265);

            auto g_z_0_xyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 266);

            auto g_z_0_xyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 267);

            auto g_z_0_xyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 268);

            auto g_z_0_xyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 269);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyzz_xxxx, g_z_0_xyyyzz_xxxy, g_z_0_xyyyzz_xxxz, g_z_0_xyyyzz_xxyy, g_z_0_xyyyzz_xxyz, g_z_0_xyyyzz_xxzz, g_z_0_xyyyzz_xyyy, g_z_0_xyyyzz_xyyz, g_z_0_xyyyzz_xyzz, g_z_0_xyyyzz_xzzz, g_z_0_xyyyzz_yyyy, g_z_0_xyyyzz_yyyz, g_z_0_xyyyzz_yyzz, g_z_0_xyyyzz_yzzz, g_z_0_xyyyzz_zzzz, g_z_0_yyyzz_xxxx, g_z_0_yyyzz_xxxxx, g_z_0_yyyzz_xxxxy, g_z_0_yyyzz_xxxxz, g_z_0_yyyzz_xxxy, g_z_0_yyyzz_xxxyy, g_z_0_yyyzz_xxxyz, g_z_0_yyyzz_xxxz, g_z_0_yyyzz_xxxzz, g_z_0_yyyzz_xxyy, g_z_0_yyyzz_xxyyy, g_z_0_yyyzz_xxyyz, g_z_0_yyyzz_xxyz, g_z_0_yyyzz_xxyzz, g_z_0_yyyzz_xxzz, g_z_0_yyyzz_xxzzz, g_z_0_yyyzz_xyyy, g_z_0_yyyzz_xyyyy, g_z_0_yyyzz_xyyyz, g_z_0_yyyzz_xyyz, g_z_0_yyyzz_xyyzz, g_z_0_yyyzz_xyzz, g_z_0_yyyzz_xyzzz, g_z_0_yyyzz_xzzz, g_z_0_yyyzz_xzzzz, g_z_0_yyyzz_yyyy, g_z_0_yyyzz_yyyz, g_z_0_yyyzz_yyzz, g_z_0_yyyzz_yzzz, g_z_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_xxxx[k] = -g_z_0_yyyzz_xxxx[k] * cd_x[k] + g_z_0_yyyzz_xxxxx[k];

                g_z_0_xyyyzz_xxxy[k] = -g_z_0_yyyzz_xxxy[k] * cd_x[k] + g_z_0_yyyzz_xxxxy[k];

                g_z_0_xyyyzz_xxxz[k] = -g_z_0_yyyzz_xxxz[k] * cd_x[k] + g_z_0_yyyzz_xxxxz[k];

                g_z_0_xyyyzz_xxyy[k] = -g_z_0_yyyzz_xxyy[k] * cd_x[k] + g_z_0_yyyzz_xxxyy[k];

                g_z_0_xyyyzz_xxyz[k] = -g_z_0_yyyzz_xxyz[k] * cd_x[k] + g_z_0_yyyzz_xxxyz[k];

                g_z_0_xyyyzz_xxzz[k] = -g_z_0_yyyzz_xxzz[k] * cd_x[k] + g_z_0_yyyzz_xxxzz[k];

                g_z_0_xyyyzz_xyyy[k] = -g_z_0_yyyzz_xyyy[k] * cd_x[k] + g_z_0_yyyzz_xxyyy[k];

                g_z_0_xyyyzz_xyyz[k] = -g_z_0_yyyzz_xyyz[k] * cd_x[k] + g_z_0_yyyzz_xxyyz[k];

                g_z_0_xyyyzz_xyzz[k] = -g_z_0_yyyzz_xyzz[k] * cd_x[k] + g_z_0_yyyzz_xxyzz[k];

                g_z_0_xyyyzz_xzzz[k] = -g_z_0_yyyzz_xzzz[k] * cd_x[k] + g_z_0_yyyzz_xxzzz[k];

                g_z_0_xyyyzz_yyyy[k] = -g_z_0_yyyzz_yyyy[k] * cd_x[k] + g_z_0_yyyzz_xyyyy[k];

                g_z_0_xyyyzz_yyyz[k] = -g_z_0_yyyzz_yyyz[k] * cd_x[k] + g_z_0_yyyzz_xyyyz[k];

                g_z_0_xyyyzz_yyzz[k] = -g_z_0_yyyzz_yyzz[k] * cd_x[k] + g_z_0_yyyzz_xyyzz[k];

                g_z_0_xyyyzz_yzzz[k] = -g_z_0_yyyzz_yzzz[k] * cd_x[k] + g_z_0_yyyzz_xyzzz[k];

                g_z_0_xyyyzz_zzzz[k] = -g_z_0_yyyzz_zzzz[k] * cd_x[k] + g_z_0_yyyzz_xzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 270);

            auto g_z_0_xyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 271);

            auto g_z_0_xyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 272);

            auto g_z_0_xyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 273);

            auto g_z_0_xyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 274);

            auto g_z_0_xyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 275);

            auto g_z_0_xyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 276);

            auto g_z_0_xyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 277);

            auto g_z_0_xyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 278);

            auto g_z_0_xyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 279);

            auto g_z_0_xyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 280);

            auto g_z_0_xyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 281);

            auto g_z_0_xyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 282);

            auto g_z_0_xyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 283);

            auto g_z_0_xyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 284);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzzz_xxxx, g_z_0_xyyzzz_xxxy, g_z_0_xyyzzz_xxxz, g_z_0_xyyzzz_xxyy, g_z_0_xyyzzz_xxyz, g_z_0_xyyzzz_xxzz, g_z_0_xyyzzz_xyyy, g_z_0_xyyzzz_xyyz, g_z_0_xyyzzz_xyzz, g_z_0_xyyzzz_xzzz, g_z_0_xyyzzz_yyyy, g_z_0_xyyzzz_yyyz, g_z_0_xyyzzz_yyzz, g_z_0_xyyzzz_yzzz, g_z_0_xyyzzz_zzzz, g_z_0_yyzzz_xxxx, g_z_0_yyzzz_xxxxx, g_z_0_yyzzz_xxxxy, g_z_0_yyzzz_xxxxz, g_z_0_yyzzz_xxxy, g_z_0_yyzzz_xxxyy, g_z_0_yyzzz_xxxyz, g_z_0_yyzzz_xxxz, g_z_0_yyzzz_xxxzz, g_z_0_yyzzz_xxyy, g_z_0_yyzzz_xxyyy, g_z_0_yyzzz_xxyyz, g_z_0_yyzzz_xxyz, g_z_0_yyzzz_xxyzz, g_z_0_yyzzz_xxzz, g_z_0_yyzzz_xxzzz, g_z_0_yyzzz_xyyy, g_z_0_yyzzz_xyyyy, g_z_0_yyzzz_xyyyz, g_z_0_yyzzz_xyyz, g_z_0_yyzzz_xyyzz, g_z_0_yyzzz_xyzz, g_z_0_yyzzz_xyzzz, g_z_0_yyzzz_xzzz, g_z_0_yyzzz_xzzzz, g_z_0_yyzzz_yyyy, g_z_0_yyzzz_yyyz, g_z_0_yyzzz_yyzz, g_z_0_yyzzz_yzzz, g_z_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_xxxx[k] = -g_z_0_yyzzz_xxxx[k] * cd_x[k] + g_z_0_yyzzz_xxxxx[k];

                g_z_0_xyyzzz_xxxy[k] = -g_z_0_yyzzz_xxxy[k] * cd_x[k] + g_z_0_yyzzz_xxxxy[k];

                g_z_0_xyyzzz_xxxz[k] = -g_z_0_yyzzz_xxxz[k] * cd_x[k] + g_z_0_yyzzz_xxxxz[k];

                g_z_0_xyyzzz_xxyy[k] = -g_z_0_yyzzz_xxyy[k] * cd_x[k] + g_z_0_yyzzz_xxxyy[k];

                g_z_0_xyyzzz_xxyz[k] = -g_z_0_yyzzz_xxyz[k] * cd_x[k] + g_z_0_yyzzz_xxxyz[k];

                g_z_0_xyyzzz_xxzz[k] = -g_z_0_yyzzz_xxzz[k] * cd_x[k] + g_z_0_yyzzz_xxxzz[k];

                g_z_0_xyyzzz_xyyy[k] = -g_z_0_yyzzz_xyyy[k] * cd_x[k] + g_z_0_yyzzz_xxyyy[k];

                g_z_0_xyyzzz_xyyz[k] = -g_z_0_yyzzz_xyyz[k] * cd_x[k] + g_z_0_yyzzz_xxyyz[k];

                g_z_0_xyyzzz_xyzz[k] = -g_z_0_yyzzz_xyzz[k] * cd_x[k] + g_z_0_yyzzz_xxyzz[k];

                g_z_0_xyyzzz_xzzz[k] = -g_z_0_yyzzz_xzzz[k] * cd_x[k] + g_z_0_yyzzz_xxzzz[k];

                g_z_0_xyyzzz_yyyy[k] = -g_z_0_yyzzz_yyyy[k] * cd_x[k] + g_z_0_yyzzz_xyyyy[k];

                g_z_0_xyyzzz_yyyz[k] = -g_z_0_yyzzz_yyyz[k] * cd_x[k] + g_z_0_yyzzz_xyyyz[k];

                g_z_0_xyyzzz_yyzz[k] = -g_z_0_yyzzz_yyzz[k] * cd_x[k] + g_z_0_yyzzz_xyyzz[k];

                g_z_0_xyyzzz_yzzz[k] = -g_z_0_yyzzz_yzzz[k] * cd_x[k] + g_z_0_yyzzz_xyzzz[k];

                g_z_0_xyyzzz_zzzz[k] = -g_z_0_yyzzz_zzzz[k] * cd_x[k] + g_z_0_yyzzz_xzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 285);

            auto g_z_0_xyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 286);

            auto g_z_0_xyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 287);

            auto g_z_0_xyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 288);

            auto g_z_0_xyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 289);

            auto g_z_0_xyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 290);

            auto g_z_0_xyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 291);

            auto g_z_0_xyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 292);

            auto g_z_0_xyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 293);

            auto g_z_0_xyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 294);

            auto g_z_0_xyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 295);

            auto g_z_0_xyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 296);

            auto g_z_0_xyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 297);

            auto g_z_0_xyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 298);

            auto g_z_0_xyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 299);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzzz_xxxx, g_z_0_xyzzzz_xxxy, g_z_0_xyzzzz_xxxz, g_z_0_xyzzzz_xxyy, g_z_0_xyzzzz_xxyz, g_z_0_xyzzzz_xxzz, g_z_0_xyzzzz_xyyy, g_z_0_xyzzzz_xyyz, g_z_0_xyzzzz_xyzz, g_z_0_xyzzzz_xzzz, g_z_0_xyzzzz_yyyy, g_z_0_xyzzzz_yyyz, g_z_0_xyzzzz_yyzz, g_z_0_xyzzzz_yzzz, g_z_0_xyzzzz_zzzz, g_z_0_yzzzz_xxxx, g_z_0_yzzzz_xxxxx, g_z_0_yzzzz_xxxxy, g_z_0_yzzzz_xxxxz, g_z_0_yzzzz_xxxy, g_z_0_yzzzz_xxxyy, g_z_0_yzzzz_xxxyz, g_z_0_yzzzz_xxxz, g_z_0_yzzzz_xxxzz, g_z_0_yzzzz_xxyy, g_z_0_yzzzz_xxyyy, g_z_0_yzzzz_xxyyz, g_z_0_yzzzz_xxyz, g_z_0_yzzzz_xxyzz, g_z_0_yzzzz_xxzz, g_z_0_yzzzz_xxzzz, g_z_0_yzzzz_xyyy, g_z_0_yzzzz_xyyyy, g_z_0_yzzzz_xyyyz, g_z_0_yzzzz_xyyz, g_z_0_yzzzz_xyyzz, g_z_0_yzzzz_xyzz, g_z_0_yzzzz_xyzzz, g_z_0_yzzzz_xzzz, g_z_0_yzzzz_xzzzz, g_z_0_yzzzz_yyyy, g_z_0_yzzzz_yyyz, g_z_0_yzzzz_yyzz, g_z_0_yzzzz_yzzz, g_z_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_xxxx[k] = -g_z_0_yzzzz_xxxx[k] * cd_x[k] + g_z_0_yzzzz_xxxxx[k];

                g_z_0_xyzzzz_xxxy[k] = -g_z_0_yzzzz_xxxy[k] * cd_x[k] + g_z_0_yzzzz_xxxxy[k];

                g_z_0_xyzzzz_xxxz[k] = -g_z_0_yzzzz_xxxz[k] * cd_x[k] + g_z_0_yzzzz_xxxxz[k];

                g_z_0_xyzzzz_xxyy[k] = -g_z_0_yzzzz_xxyy[k] * cd_x[k] + g_z_0_yzzzz_xxxyy[k];

                g_z_0_xyzzzz_xxyz[k] = -g_z_0_yzzzz_xxyz[k] * cd_x[k] + g_z_0_yzzzz_xxxyz[k];

                g_z_0_xyzzzz_xxzz[k] = -g_z_0_yzzzz_xxzz[k] * cd_x[k] + g_z_0_yzzzz_xxxzz[k];

                g_z_0_xyzzzz_xyyy[k] = -g_z_0_yzzzz_xyyy[k] * cd_x[k] + g_z_0_yzzzz_xxyyy[k];

                g_z_0_xyzzzz_xyyz[k] = -g_z_0_yzzzz_xyyz[k] * cd_x[k] + g_z_0_yzzzz_xxyyz[k];

                g_z_0_xyzzzz_xyzz[k] = -g_z_0_yzzzz_xyzz[k] * cd_x[k] + g_z_0_yzzzz_xxyzz[k];

                g_z_0_xyzzzz_xzzz[k] = -g_z_0_yzzzz_xzzz[k] * cd_x[k] + g_z_0_yzzzz_xxzzz[k];

                g_z_0_xyzzzz_yyyy[k] = -g_z_0_yzzzz_yyyy[k] * cd_x[k] + g_z_0_yzzzz_xyyyy[k];

                g_z_0_xyzzzz_yyyz[k] = -g_z_0_yzzzz_yyyz[k] * cd_x[k] + g_z_0_yzzzz_xyyyz[k];

                g_z_0_xyzzzz_yyzz[k] = -g_z_0_yzzzz_yyzz[k] * cd_x[k] + g_z_0_yzzzz_xyyzz[k];

                g_z_0_xyzzzz_yzzz[k] = -g_z_0_yzzzz_yzzz[k] * cd_x[k] + g_z_0_yzzzz_xyzzz[k];

                g_z_0_xyzzzz_zzzz[k] = -g_z_0_yzzzz_zzzz[k] * cd_x[k] + g_z_0_yzzzz_xzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 300);

            auto g_z_0_xzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 301);

            auto g_z_0_xzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 302);

            auto g_z_0_xzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 303);

            auto g_z_0_xzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 304);

            auto g_z_0_xzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 305);

            auto g_z_0_xzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 306);

            auto g_z_0_xzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 307);

            auto g_z_0_xzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 308);

            auto g_z_0_xzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 309);

            auto g_z_0_xzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 310);

            auto g_z_0_xzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 311);

            auto g_z_0_xzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 312);

            auto g_z_0_xzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 313);

            auto g_z_0_xzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzzz_xxxx, g_z_0_xzzzzz_xxxy, g_z_0_xzzzzz_xxxz, g_z_0_xzzzzz_xxyy, g_z_0_xzzzzz_xxyz, g_z_0_xzzzzz_xxzz, g_z_0_xzzzzz_xyyy, g_z_0_xzzzzz_xyyz, g_z_0_xzzzzz_xyzz, g_z_0_xzzzzz_xzzz, g_z_0_xzzzzz_yyyy, g_z_0_xzzzzz_yyyz, g_z_0_xzzzzz_yyzz, g_z_0_xzzzzz_yzzz, g_z_0_xzzzzz_zzzz, g_z_0_zzzzz_xxxx, g_z_0_zzzzz_xxxxx, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_yyyy, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_xxxx[k] = -g_z_0_zzzzz_xxxx[k] * cd_x[k] + g_z_0_zzzzz_xxxxx[k];

                g_z_0_xzzzzz_xxxy[k] = -g_z_0_zzzzz_xxxy[k] * cd_x[k] + g_z_0_zzzzz_xxxxy[k];

                g_z_0_xzzzzz_xxxz[k] = -g_z_0_zzzzz_xxxz[k] * cd_x[k] + g_z_0_zzzzz_xxxxz[k];

                g_z_0_xzzzzz_xxyy[k] = -g_z_0_zzzzz_xxyy[k] * cd_x[k] + g_z_0_zzzzz_xxxyy[k];

                g_z_0_xzzzzz_xxyz[k] = -g_z_0_zzzzz_xxyz[k] * cd_x[k] + g_z_0_zzzzz_xxxyz[k];

                g_z_0_xzzzzz_xxzz[k] = -g_z_0_zzzzz_xxzz[k] * cd_x[k] + g_z_0_zzzzz_xxxzz[k];

                g_z_0_xzzzzz_xyyy[k] = -g_z_0_zzzzz_xyyy[k] * cd_x[k] + g_z_0_zzzzz_xxyyy[k];

                g_z_0_xzzzzz_xyyz[k] = -g_z_0_zzzzz_xyyz[k] * cd_x[k] + g_z_0_zzzzz_xxyyz[k];

                g_z_0_xzzzzz_xyzz[k] = -g_z_0_zzzzz_xyzz[k] * cd_x[k] + g_z_0_zzzzz_xxyzz[k];

                g_z_0_xzzzzz_xzzz[k] = -g_z_0_zzzzz_xzzz[k] * cd_x[k] + g_z_0_zzzzz_xxzzz[k];

                g_z_0_xzzzzz_yyyy[k] = -g_z_0_zzzzz_yyyy[k] * cd_x[k] + g_z_0_zzzzz_xyyyy[k];

                g_z_0_xzzzzz_yyyz[k] = -g_z_0_zzzzz_yyyz[k] * cd_x[k] + g_z_0_zzzzz_xyyyz[k];

                g_z_0_xzzzzz_yyzz[k] = -g_z_0_zzzzz_yyzz[k] * cd_x[k] + g_z_0_zzzzz_xyyzz[k];

                g_z_0_xzzzzz_yzzz[k] = -g_z_0_zzzzz_yzzz[k] * cd_x[k] + g_z_0_zzzzz_xyzzz[k];

                g_z_0_xzzzzz_zzzz[k] = -g_z_0_zzzzz_zzzz[k] * cd_x[k] + g_z_0_zzzzz_xzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 315);

            auto g_z_0_yyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 316);

            auto g_z_0_yyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 317);

            auto g_z_0_yyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 318);

            auto g_z_0_yyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 319);

            auto g_z_0_yyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 320);

            auto g_z_0_yyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 321);

            auto g_z_0_yyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 322);

            auto g_z_0_yyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 323);

            auto g_z_0_yyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 324);

            auto g_z_0_yyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 325);

            auto g_z_0_yyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 326);

            auto g_z_0_yyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 327);

            auto g_z_0_yyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 328);

            auto g_z_0_yyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 329);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyy_xxxx, g_z_0_yyyyy_xxxxy, g_z_0_yyyyy_xxxy, g_z_0_yyyyy_xxxyy, g_z_0_yyyyy_xxxyz, g_z_0_yyyyy_xxxz, g_z_0_yyyyy_xxyy, g_z_0_yyyyy_xxyyy, g_z_0_yyyyy_xxyyz, g_z_0_yyyyy_xxyz, g_z_0_yyyyy_xxyzz, g_z_0_yyyyy_xxzz, g_z_0_yyyyy_xyyy, g_z_0_yyyyy_xyyyy, g_z_0_yyyyy_xyyyz, g_z_0_yyyyy_xyyz, g_z_0_yyyyy_xyyzz, g_z_0_yyyyy_xyzz, g_z_0_yyyyy_xyzzz, g_z_0_yyyyy_xzzz, g_z_0_yyyyy_yyyy, g_z_0_yyyyy_yyyyy, g_z_0_yyyyy_yyyyz, g_z_0_yyyyy_yyyz, g_z_0_yyyyy_yyyzz, g_z_0_yyyyy_yyzz, g_z_0_yyyyy_yyzzz, g_z_0_yyyyy_yzzz, g_z_0_yyyyy_yzzzz, g_z_0_yyyyy_zzzz, g_z_0_yyyyyy_xxxx, g_z_0_yyyyyy_xxxy, g_z_0_yyyyyy_xxxz, g_z_0_yyyyyy_xxyy, g_z_0_yyyyyy_xxyz, g_z_0_yyyyyy_xxzz, g_z_0_yyyyyy_xyyy, g_z_0_yyyyyy_xyyz, g_z_0_yyyyyy_xyzz, g_z_0_yyyyyy_xzzz, g_z_0_yyyyyy_yyyy, g_z_0_yyyyyy_yyyz, g_z_0_yyyyyy_yyzz, g_z_0_yyyyyy_yzzz, g_z_0_yyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_xxxx[k] = -g_z_0_yyyyy_xxxx[k] * cd_y[k] + g_z_0_yyyyy_xxxxy[k];

                g_z_0_yyyyyy_xxxy[k] = -g_z_0_yyyyy_xxxy[k] * cd_y[k] + g_z_0_yyyyy_xxxyy[k];

                g_z_0_yyyyyy_xxxz[k] = -g_z_0_yyyyy_xxxz[k] * cd_y[k] + g_z_0_yyyyy_xxxyz[k];

                g_z_0_yyyyyy_xxyy[k] = -g_z_0_yyyyy_xxyy[k] * cd_y[k] + g_z_0_yyyyy_xxyyy[k];

                g_z_0_yyyyyy_xxyz[k] = -g_z_0_yyyyy_xxyz[k] * cd_y[k] + g_z_0_yyyyy_xxyyz[k];

                g_z_0_yyyyyy_xxzz[k] = -g_z_0_yyyyy_xxzz[k] * cd_y[k] + g_z_0_yyyyy_xxyzz[k];

                g_z_0_yyyyyy_xyyy[k] = -g_z_0_yyyyy_xyyy[k] * cd_y[k] + g_z_0_yyyyy_xyyyy[k];

                g_z_0_yyyyyy_xyyz[k] = -g_z_0_yyyyy_xyyz[k] * cd_y[k] + g_z_0_yyyyy_xyyyz[k];

                g_z_0_yyyyyy_xyzz[k] = -g_z_0_yyyyy_xyzz[k] * cd_y[k] + g_z_0_yyyyy_xyyzz[k];

                g_z_0_yyyyyy_xzzz[k] = -g_z_0_yyyyy_xzzz[k] * cd_y[k] + g_z_0_yyyyy_xyzzz[k];

                g_z_0_yyyyyy_yyyy[k] = -g_z_0_yyyyy_yyyy[k] * cd_y[k] + g_z_0_yyyyy_yyyyy[k];

                g_z_0_yyyyyy_yyyz[k] = -g_z_0_yyyyy_yyyz[k] * cd_y[k] + g_z_0_yyyyy_yyyyz[k];

                g_z_0_yyyyyy_yyzz[k] = -g_z_0_yyyyy_yyzz[k] * cd_y[k] + g_z_0_yyyyy_yyyzz[k];

                g_z_0_yyyyyy_yzzz[k] = -g_z_0_yyyyy_yzzz[k] * cd_y[k] + g_z_0_yyyyy_yyzzz[k];

                g_z_0_yyyyyy_zzzz[k] = -g_z_0_yyyyy_zzzz[k] * cd_y[k] + g_z_0_yyyyy_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 330);

            auto g_z_0_yyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 331);

            auto g_z_0_yyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 332);

            auto g_z_0_yyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 333);

            auto g_z_0_yyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 334);

            auto g_z_0_yyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 335);

            auto g_z_0_yyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 336);

            auto g_z_0_yyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 337);

            auto g_z_0_yyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 338);

            auto g_z_0_yyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 339);

            auto g_z_0_yyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 340);

            auto g_z_0_yyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 341);

            auto g_z_0_yyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 342);

            auto g_z_0_yyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 343);

            auto g_z_0_yyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 344);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyyz_xxxx, g_z_0_yyyyyz_xxxy, g_z_0_yyyyyz_xxxz, g_z_0_yyyyyz_xxyy, g_z_0_yyyyyz_xxyz, g_z_0_yyyyyz_xxzz, g_z_0_yyyyyz_xyyy, g_z_0_yyyyyz_xyyz, g_z_0_yyyyyz_xyzz, g_z_0_yyyyyz_xzzz, g_z_0_yyyyyz_yyyy, g_z_0_yyyyyz_yyyz, g_z_0_yyyyyz_yyzz, g_z_0_yyyyyz_yzzz, g_z_0_yyyyyz_zzzz, g_z_0_yyyyz_xxxx, g_z_0_yyyyz_xxxxy, g_z_0_yyyyz_xxxy, g_z_0_yyyyz_xxxyy, g_z_0_yyyyz_xxxyz, g_z_0_yyyyz_xxxz, g_z_0_yyyyz_xxyy, g_z_0_yyyyz_xxyyy, g_z_0_yyyyz_xxyyz, g_z_0_yyyyz_xxyz, g_z_0_yyyyz_xxyzz, g_z_0_yyyyz_xxzz, g_z_0_yyyyz_xyyy, g_z_0_yyyyz_xyyyy, g_z_0_yyyyz_xyyyz, g_z_0_yyyyz_xyyz, g_z_0_yyyyz_xyyzz, g_z_0_yyyyz_xyzz, g_z_0_yyyyz_xyzzz, g_z_0_yyyyz_xzzz, g_z_0_yyyyz_yyyy, g_z_0_yyyyz_yyyyy, g_z_0_yyyyz_yyyyz, g_z_0_yyyyz_yyyz, g_z_0_yyyyz_yyyzz, g_z_0_yyyyz_yyzz, g_z_0_yyyyz_yyzzz, g_z_0_yyyyz_yzzz, g_z_0_yyyyz_yzzzz, g_z_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_xxxx[k] = -g_z_0_yyyyz_xxxx[k] * cd_y[k] + g_z_0_yyyyz_xxxxy[k];

                g_z_0_yyyyyz_xxxy[k] = -g_z_0_yyyyz_xxxy[k] * cd_y[k] + g_z_0_yyyyz_xxxyy[k];

                g_z_0_yyyyyz_xxxz[k] = -g_z_0_yyyyz_xxxz[k] * cd_y[k] + g_z_0_yyyyz_xxxyz[k];

                g_z_0_yyyyyz_xxyy[k] = -g_z_0_yyyyz_xxyy[k] * cd_y[k] + g_z_0_yyyyz_xxyyy[k];

                g_z_0_yyyyyz_xxyz[k] = -g_z_0_yyyyz_xxyz[k] * cd_y[k] + g_z_0_yyyyz_xxyyz[k];

                g_z_0_yyyyyz_xxzz[k] = -g_z_0_yyyyz_xxzz[k] * cd_y[k] + g_z_0_yyyyz_xxyzz[k];

                g_z_0_yyyyyz_xyyy[k] = -g_z_0_yyyyz_xyyy[k] * cd_y[k] + g_z_0_yyyyz_xyyyy[k];

                g_z_0_yyyyyz_xyyz[k] = -g_z_0_yyyyz_xyyz[k] * cd_y[k] + g_z_0_yyyyz_xyyyz[k];

                g_z_0_yyyyyz_xyzz[k] = -g_z_0_yyyyz_xyzz[k] * cd_y[k] + g_z_0_yyyyz_xyyzz[k];

                g_z_0_yyyyyz_xzzz[k] = -g_z_0_yyyyz_xzzz[k] * cd_y[k] + g_z_0_yyyyz_xyzzz[k];

                g_z_0_yyyyyz_yyyy[k] = -g_z_0_yyyyz_yyyy[k] * cd_y[k] + g_z_0_yyyyz_yyyyy[k];

                g_z_0_yyyyyz_yyyz[k] = -g_z_0_yyyyz_yyyz[k] * cd_y[k] + g_z_0_yyyyz_yyyyz[k];

                g_z_0_yyyyyz_yyzz[k] = -g_z_0_yyyyz_yyzz[k] * cd_y[k] + g_z_0_yyyyz_yyyzz[k];

                g_z_0_yyyyyz_yzzz[k] = -g_z_0_yyyyz_yzzz[k] * cd_y[k] + g_z_0_yyyyz_yyzzz[k];

                g_z_0_yyyyyz_zzzz[k] = -g_z_0_yyyyz_zzzz[k] * cd_y[k] + g_z_0_yyyyz_yzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 345);

            auto g_z_0_yyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 346);

            auto g_z_0_yyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 347);

            auto g_z_0_yyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 348);

            auto g_z_0_yyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 349);

            auto g_z_0_yyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 350);

            auto g_z_0_yyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 351);

            auto g_z_0_yyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 352);

            auto g_z_0_yyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 353);

            auto g_z_0_yyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 354);

            auto g_z_0_yyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 355);

            auto g_z_0_yyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 356);

            auto g_z_0_yyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 357);

            auto g_z_0_yyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 358);

            auto g_z_0_yyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 359);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyzz_xxxx, g_z_0_yyyyzz_xxxy, g_z_0_yyyyzz_xxxz, g_z_0_yyyyzz_xxyy, g_z_0_yyyyzz_xxyz, g_z_0_yyyyzz_xxzz, g_z_0_yyyyzz_xyyy, g_z_0_yyyyzz_xyyz, g_z_0_yyyyzz_xyzz, g_z_0_yyyyzz_xzzz, g_z_0_yyyyzz_yyyy, g_z_0_yyyyzz_yyyz, g_z_0_yyyyzz_yyzz, g_z_0_yyyyzz_yzzz, g_z_0_yyyyzz_zzzz, g_z_0_yyyzz_xxxx, g_z_0_yyyzz_xxxxy, g_z_0_yyyzz_xxxy, g_z_0_yyyzz_xxxyy, g_z_0_yyyzz_xxxyz, g_z_0_yyyzz_xxxz, g_z_0_yyyzz_xxyy, g_z_0_yyyzz_xxyyy, g_z_0_yyyzz_xxyyz, g_z_0_yyyzz_xxyz, g_z_0_yyyzz_xxyzz, g_z_0_yyyzz_xxzz, g_z_0_yyyzz_xyyy, g_z_0_yyyzz_xyyyy, g_z_0_yyyzz_xyyyz, g_z_0_yyyzz_xyyz, g_z_0_yyyzz_xyyzz, g_z_0_yyyzz_xyzz, g_z_0_yyyzz_xyzzz, g_z_0_yyyzz_xzzz, g_z_0_yyyzz_yyyy, g_z_0_yyyzz_yyyyy, g_z_0_yyyzz_yyyyz, g_z_0_yyyzz_yyyz, g_z_0_yyyzz_yyyzz, g_z_0_yyyzz_yyzz, g_z_0_yyyzz_yyzzz, g_z_0_yyyzz_yzzz, g_z_0_yyyzz_yzzzz, g_z_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_xxxx[k] = -g_z_0_yyyzz_xxxx[k] * cd_y[k] + g_z_0_yyyzz_xxxxy[k];

                g_z_0_yyyyzz_xxxy[k] = -g_z_0_yyyzz_xxxy[k] * cd_y[k] + g_z_0_yyyzz_xxxyy[k];

                g_z_0_yyyyzz_xxxz[k] = -g_z_0_yyyzz_xxxz[k] * cd_y[k] + g_z_0_yyyzz_xxxyz[k];

                g_z_0_yyyyzz_xxyy[k] = -g_z_0_yyyzz_xxyy[k] * cd_y[k] + g_z_0_yyyzz_xxyyy[k];

                g_z_0_yyyyzz_xxyz[k] = -g_z_0_yyyzz_xxyz[k] * cd_y[k] + g_z_0_yyyzz_xxyyz[k];

                g_z_0_yyyyzz_xxzz[k] = -g_z_0_yyyzz_xxzz[k] * cd_y[k] + g_z_0_yyyzz_xxyzz[k];

                g_z_0_yyyyzz_xyyy[k] = -g_z_0_yyyzz_xyyy[k] * cd_y[k] + g_z_0_yyyzz_xyyyy[k];

                g_z_0_yyyyzz_xyyz[k] = -g_z_0_yyyzz_xyyz[k] * cd_y[k] + g_z_0_yyyzz_xyyyz[k];

                g_z_0_yyyyzz_xyzz[k] = -g_z_0_yyyzz_xyzz[k] * cd_y[k] + g_z_0_yyyzz_xyyzz[k];

                g_z_0_yyyyzz_xzzz[k] = -g_z_0_yyyzz_xzzz[k] * cd_y[k] + g_z_0_yyyzz_xyzzz[k];

                g_z_0_yyyyzz_yyyy[k] = -g_z_0_yyyzz_yyyy[k] * cd_y[k] + g_z_0_yyyzz_yyyyy[k];

                g_z_0_yyyyzz_yyyz[k] = -g_z_0_yyyzz_yyyz[k] * cd_y[k] + g_z_0_yyyzz_yyyyz[k];

                g_z_0_yyyyzz_yyzz[k] = -g_z_0_yyyzz_yyzz[k] * cd_y[k] + g_z_0_yyyzz_yyyzz[k];

                g_z_0_yyyyzz_yzzz[k] = -g_z_0_yyyzz_yzzz[k] * cd_y[k] + g_z_0_yyyzz_yyzzz[k];

                g_z_0_yyyyzz_zzzz[k] = -g_z_0_yyyzz_zzzz[k] * cd_y[k] + g_z_0_yyyzz_yzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 360);

            auto g_z_0_yyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 361);

            auto g_z_0_yyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 362);

            auto g_z_0_yyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 363);

            auto g_z_0_yyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 364);

            auto g_z_0_yyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 365);

            auto g_z_0_yyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 366);

            auto g_z_0_yyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 367);

            auto g_z_0_yyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 368);

            auto g_z_0_yyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 369);

            auto g_z_0_yyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 370);

            auto g_z_0_yyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 371);

            auto g_z_0_yyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 372);

            auto g_z_0_yyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 373);

            auto g_z_0_yyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 374);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzzz_xxxx, g_z_0_yyyzzz_xxxy, g_z_0_yyyzzz_xxxz, g_z_0_yyyzzz_xxyy, g_z_0_yyyzzz_xxyz, g_z_0_yyyzzz_xxzz, g_z_0_yyyzzz_xyyy, g_z_0_yyyzzz_xyyz, g_z_0_yyyzzz_xyzz, g_z_0_yyyzzz_xzzz, g_z_0_yyyzzz_yyyy, g_z_0_yyyzzz_yyyz, g_z_0_yyyzzz_yyzz, g_z_0_yyyzzz_yzzz, g_z_0_yyyzzz_zzzz, g_z_0_yyzzz_xxxx, g_z_0_yyzzz_xxxxy, g_z_0_yyzzz_xxxy, g_z_0_yyzzz_xxxyy, g_z_0_yyzzz_xxxyz, g_z_0_yyzzz_xxxz, g_z_0_yyzzz_xxyy, g_z_0_yyzzz_xxyyy, g_z_0_yyzzz_xxyyz, g_z_0_yyzzz_xxyz, g_z_0_yyzzz_xxyzz, g_z_0_yyzzz_xxzz, g_z_0_yyzzz_xyyy, g_z_0_yyzzz_xyyyy, g_z_0_yyzzz_xyyyz, g_z_0_yyzzz_xyyz, g_z_0_yyzzz_xyyzz, g_z_0_yyzzz_xyzz, g_z_0_yyzzz_xyzzz, g_z_0_yyzzz_xzzz, g_z_0_yyzzz_yyyy, g_z_0_yyzzz_yyyyy, g_z_0_yyzzz_yyyyz, g_z_0_yyzzz_yyyz, g_z_0_yyzzz_yyyzz, g_z_0_yyzzz_yyzz, g_z_0_yyzzz_yyzzz, g_z_0_yyzzz_yzzz, g_z_0_yyzzz_yzzzz, g_z_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_xxxx[k] = -g_z_0_yyzzz_xxxx[k] * cd_y[k] + g_z_0_yyzzz_xxxxy[k];

                g_z_0_yyyzzz_xxxy[k] = -g_z_0_yyzzz_xxxy[k] * cd_y[k] + g_z_0_yyzzz_xxxyy[k];

                g_z_0_yyyzzz_xxxz[k] = -g_z_0_yyzzz_xxxz[k] * cd_y[k] + g_z_0_yyzzz_xxxyz[k];

                g_z_0_yyyzzz_xxyy[k] = -g_z_0_yyzzz_xxyy[k] * cd_y[k] + g_z_0_yyzzz_xxyyy[k];

                g_z_0_yyyzzz_xxyz[k] = -g_z_0_yyzzz_xxyz[k] * cd_y[k] + g_z_0_yyzzz_xxyyz[k];

                g_z_0_yyyzzz_xxzz[k] = -g_z_0_yyzzz_xxzz[k] * cd_y[k] + g_z_0_yyzzz_xxyzz[k];

                g_z_0_yyyzzz_xyyy[k] = -g_z_0_yyzzz_xyyy[k] * cd_y[k] + g_z_0_yyzzz_xyyyy[k];

                g_z_0_yyyzzz_xyyz[k] = -g_z_0_yyzzz_xyyz[k] * cd_y[k] + g_z_0_yyzzz_xyyyz[k];

                g_z_0_yyyzzz_xyzz[k] = -g_z_0_yyzzz_xyzz[k] * cd_y[k] + g_z_0_yyzzz_xyyzz[k];

                g_z_0_yyyzzz_xzzz[k] = -g_z_0_yyzzz_xzzz[k] * cd_y[k] + g_z_0_yyzzz_xyzzz[k];

                g_z_0_yyyzzz_yyyy[k] = -g_z_0_yyzzz_yyyy[k] * cd_y[k] + g_z_0_yyzzz_yyyyy[k];

                g_z_0_yyyzzz_yyyz[k] = -g_z_0_yyzzz_yyyz[k] * cd_y[k] + g_z_0_yyzzz_yyyyz[k];

                g_z_0_yyyzzz_yyzz[k] = -g_z_0_yyzzz_yyzz[k] * cd_y[k] + g_z_0_yyzzz_yyyzz[k];

                g_z_0_yyyzzz_yzzz[k] = -g_z_0_yyzzz_yzzz[k] * cd_y[k] + g_z_0_yyzzz_yyzzz[k];

                g_z_0_yyyzzz_zzzz[k] = -g_z_0_yyzzz_zzzz[k] * cd_y[k] + g_z_0_yyzzz_yzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 375);

            auto g_z_0_yyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 376);

            auto g_z_0_yyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 377);

            auto g_z_0_yyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 378);

            auto g_z_0_yyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 379);

            auto g_z_0_yyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 380);

            auto g_z_0_yyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 381);

            auto g_z_0_yyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 382);

            auto g_z_0_yyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 383);

            auto g_z_0_yyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 384);

            auto g_z_0_yyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 385);

            auto g_z_0_yyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 386);

            auto g_z_0_yyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 387);

            auto g_z_0_yyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 388);

            auto g_z_0_yyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 389);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzzz_xxxx, g_z_0_yyzzzz_xxxy, g_z_0_yyzzzz_xxxz, g_z_0_yyzzzz_xxyy, g_z_0_yyzzzz_xxyz, g_z_0_yyzzzz_xxzz, g_z_0_yyzzzz_xyyy, g_z_0_yyzzzz_xyyz, g_z_0_yyzzzz_xyzz, g_z_0_yyzzzz_xzzz, g_z_0_yyzzzz_yyyy, g_z_0_yyzzzz_yyyz, g_z_0_yyzzzz_yyzz, g_z_0_yyzzzz_yzzz, g_z_0_yyzzzz_zzzz, g_z_0_yzzzz_xxxx, g_z_0_yzzzz_xxxxy, g_z_0_yzzzz_xxxy, g_z_0_yzzzz_xxxyy, g_z_0_yzzzz_xxxyz, g_z_0_yzzzz_xxxz, g_z_0_yzzzz_xxyy, g_z_0_yzzzz_xxyyy, g_z_0_yzzzz_xxyyz, g_z_0_yzzzz_xxyz, g_z_0_yzzzz_xxyzz, g_z_0_yzzzz_xxzz, g_z_0_yzzzz_xyyy, g_z_0_yzzzz_xyyyy, g_z_0_yzzzz_xyyyz, g_z_0_yzzzz_xyyz, g_z_0_yzzzz_xyyzz, g_z_0_yzzzz_xyzz, g_z_0_yzzzz_xyzzz, g_z_0_yzzzz_xzzz, g_z_0_yzzzz_yyyy, g_z_0_yzzzz_yyyyy, g_z_0_yzzzz_yyyyz, g_z_0_yzzzz_yyyz, g_z_0_yzzzz_yyyzz, g_z_0_yzzzz_yyzz, g_z_0_yzzzz_yyzzz, g_z_0_yzzzz_yzzz, g_z_0_yzzzz_yzzzz, g_z_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_xxxx[k] = -g_z_0_yzzzz_xxxx[k] * cd_y[k] + g_z_0_yzzzz_xxxxy[k];

                g_z_0_yyzzzz_xxxy[k] = -g_z_0_yzzzz_xxxy[k] * cd_y[k] + g_z_0_yzzzz_xxxyy[k];

                g_z_0_yyzzzz_xxxz[k] = -g_z_0_yzzzz_xxxz[k] * cd_y[k] + g_z_0_yzzzz_xxxyz[k];

                g_z_0_yyzzzz_xxyy[k] = -g_z_0_yzzzz_xxyy[k] * cd_y[k] + g_z_0_yzzzz_xxyyy[k];

                g_z_0_yyzzzz_xxyz[k] = -g_z_0_yzzzz_xxyz[k] * cd_y[k] + g_z_0_yzzzz_xxyyz[k];

                g_z_0_yyzzzz_xxzz[k] = -g_z_0_yzzzz_xxzz[k] * cd_y[k] + g_z_0_yzzzz_xxyzz[k];

                g_z_0_yyzzzz_xyyy[k] = -g_z_0_yzzzz_xyyy[k] * cd_y[k] + g_z_0_yzzzz_xyyyy[k];

                g_z_0_yyzzzz_xyyz[k] = -g_z_0_yzzzz_xyyz[k] * cd_y[k] + g_z_0_yzzzz_xyyyz[k];

                g_z_0_yyzzzz_xyzz[k] = -g_z_0_yzzzz_xyzz[k] * cd_y[k] + g_z_0_yzzzz_xyyzz[k];

                g_z_0_yyzzzz_xzzz[k] = -g_z_0_yzzzz_xzzz[k] * cd_y[k] + g_z_0_yzzzz_xyzzz[k];

                g_z_0_yyzzzz_yyyy[k] = -g_z_0_yzzzz_yyyy[k] * cd_y[k] + g_z_0_yzzzz_yyyyy[k];

                g_z_0_yyzzzz_yyyz[k] = -g_z_0_yzzzz_yyyz[k] * cd_y[k] + g_z_0_yzzzz_yyyyz[k];

                g_z_0_yyzzzz_yyzz[k] = -g_z_0_yzzzz_yyzz[k] * cd_y[k] + g_z_0_yzzzz_yyyzz[k];

                g_z_0_yyzzzz_yzzz[k] = -g_z_0_yzzzz_yzzz[k] * cd_y[k] + g_z_0_yzzzz_yyzzz[k];

                g_z_0_yyzzzz_zzzz[k] = -g_z_0_yzzzz_zzzz[k] * cd_y[k] + g_z_0_yzzzz_yzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 390);

            auto g_z_0_yzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 391);

            auto g_z_0_yzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 392);

            auto g_z_0_yzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 393);

            auto g_z_0_yzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 394);

            auto g_z_0_yzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 395);

            auto g_z_0_yzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 396);

            auto g_z_0_yzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 397);

            auto g_z_0_yzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 398);

            auto g_z_0_yzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 399);

            auto g_z_0_yzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 400);

            auto g_z_0_yzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 401);

            auto g_z_0_yzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 402);

            auto g_z_0_yzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 403);

            auto g_z_0_yzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 404);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzzz_xxxx, g_z_0_yzzzzz_xxxy, g_z_0_yzzzzz_xxxz, g_z_0_yzzzzz_xxyy, g_z_0_yzzzzz_xxyz, g_z_0_yzzzzz_xxzz, g_z_0_yzzzzz_xyyy, g_z_0_yzzzzz_xyyz, g_z_0_yzzzzz_xyzz, g_z_0_yzzzzz_xzzz, g_z_0_yzzzzz_yyyy, g_z_0_yzzzzz_yyyz, g_z_0_yzzzzz_yyzz, g_z_0_yzzzzz_yzzz, g_z_0_yzzzzz_zzzz, g_z_0_zzzzz_xxxx, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_yyyy, g_z_0_zzzzz_yyyyy, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_xxxx[k] = -g_z_0_zzzzz_xxxx[k] * cd_y[k] + g_z_0_zzzzz_xxxxy[k];

                g_z_0_yzzzzz_xxxy[k] = -g_z_0_zzzzz_xxxy[k] * cd_y[k] + g_z_0_zzzzz_xxxyy[k];

                g_z_0_yzzzzz_xxxz[k] = -g_z_0_zzzzz_xxxz[k] * cd_y[k] + g_z_0_zzzzz_xxxyz[k];

                g_z_0_yzzzzz_xxyy[k] = -g_z_0_zzzzz_xxyy[k] * cd_y[k] + g_z_0_zzzzz_xxyyy[k];

                g_z_0_yzzzzz_xxyz[k] = -g_z_0_zzzzz_xxyz[k] * cd_y[k] + g_z_0_zzzzz_xxyyz[k];

                g_z_0_yzzzzz_xxzz[k] = -g_z_0_zzzzz_xxzz[k] * cd_y[k] + g_z_0_zzzzz_xxyzz[k];

                g_z_0_yzzzzz_xyyy[k] = -g_z_0_zzzzz_xyyy[k] * cd_y[k] + g_z_0_zzzzz_xyyyy[k];

                g_z_0_yzzzzz_xyyz[k] = -g_z_0_zzzzz_xyyz[k] * cd_y[k] + g_z_0_zzzzz_xyyyz[k];

                g_z_0_yzzzzz_xyzz[k] = -g_z_0_zzzzz_xyzz[k] * cd_y[k] + g_z_0_zzzzz_xyyzz[k];

                g_z_0_yzzzzz_xzzz[k] = -g_z_0_zzzzz_xzzz[k] * cd_y[k] + g_z_0_zzzzz_xyzzz[k];

                g_z_0_yzzzzz_yyyy[k] = -g_z_0_zzzzz_yyyy[k] * cd_y[k] + g_z_0_zzzzz_yyyyy[k];

                g_z_0_yzzzzz_yyyz[k] = -g_z_0_zzzzz_yyyz[k] * cd_y[k] + g_z_0_zzzzz_yyyyz[k];

                g_z_0_yzzzzz_yyzz[k] = -g_z_0_zzzzz_yyzz[k] * cd_y[k] + g_z_0_zzzzz_yyyzz[k];

                g_z_0_yzzzzz_yzzz[k] = -g_z_0_zzzzz_yzzz[k] * cd_y[k] + g_z_0_zzzzz_yyzzz[k];

                g_z_0_yzzzzz_zzzz[k] = -g_z_0_zzzzz_zzzz[k] * cd_y[k] + g_z_0_zzzzz_yzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 405);

            auto g_z_0_zzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 406);

            auto g_z_0_zzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 407);

            auto g_z_0_zzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 408);

            auto g_z_0_zzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 409);

            auto g_z_0_zzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 410);

            auto g_z_0_zzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 411);

            auto g_z_0_zzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 412);

            auto g_z_0_zzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 413);

            auto g_z_0_zzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 414);

            auto g_z_0_zzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 415);

            auto g_z_0_zzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 416);

            auto g_z_0_zzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 417);

            auto g_z_0_zzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 418);

            auto g_z_0_zzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 840 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_z, g_z_0_zzzzz_xxxx, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_yyyy, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_zzzz, g_z_0_zzzzz_zzzzz, g_z_0_zzzzzz_xxxx, g_z_0_zzzzzz_xxxy, g_z_0_zzzzzz_xxxz, g_z_0_zzzzzz_xxyy, g_z_0_zzzzzz_xxyz, g_z_0_zzzzzz_xxzz, g_z_0_zzzzzz_xyyy, g_z_0_zzzzzz_xyyz, g_z_0_zzzzzz_xyzz, g_z_0_zzzzzz_xzzz, g_z_0_zzzzzz_yyyy, g_z_0_zzzzzz_yyyz, g_z_0_zzzzzz_yyzz, g_z_0_zzzzzz_yzzz, g_z_0_zzzzzz_zzzz, g_zzzzz_xxxx, g_zzzzz_xxxy, g_zzzzz_xxxz, g_zzzzz_xxyy, g_zzzzz_xxyz, g_zzzzz_xxzz, g_zzzzz_xyyy, g_zzzzz_xyyz, g_zzzzz_xyzz, g_zzzzz_xzzz, g_zzzzz_yyyy, g_zzzzz_yyyz, g_zzzzz_yyzz, g_zzzzz_yzzz, g_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_xxxx[k] = -g_zzzzz_xxxx[k] - g_z_0_zzzzz_xxxx[k] * cd_z[k] + g_z_0_zzzzz_xxxxz[k];

                g_z_0_zzzzzz_xxxy[k] = -g_zzzzz_xxxy[k] - g_z_0_zzzzz_xxxy[k] * cd_z[k] + g_z_0_zzzzz_xxxyz[k];

                g_z_0_zzzzzz_xxxz[k] = -g_zzzzz_xxxz[k] - g_z_0_zzzzz_xxxz[k] * cd_z[k] + g_z_0_zzzzz_xxxzz[k];

                g_z_0_zzzzzz_xxyy[k] = -g_zzzzz_xxyy[k] - g_z_0_zzzzz_xxyy[k] * cd_z[k] + g_z_0_zzzzz_xxyyz[k];

                g_z_0_zzzzzz_xxyz[k] = -g_zzzzz_xxyz[k] - g_z_0_zzzzz_xxyz[k] * cd_z[k] + g_z_0_zzzzz_xxyzz[k];

                g_z_0_zzzzzz_xxzz[k] = -g_zzzzz_xxzz[k] - g_z_0_zzzzz_xxzz[k] * cd_z[k] + g_z_0_zzzzz_xxzzz[k];

                g_z_0_zzzzzz_xyyy[k] = -g_zzzzz_xyyy[k] - g_z_0_zzzzz_xyyy[k] * cd_z[k] + g_z_0_zzzzz_xyyyz[k];

                g_z_0_zzzzzz_xyyz[k] = -g_zzzzz_xyyz[k] - g_z_0_zzzzz_xyyz[k] * cd_z[k] + g_z_0_zzzzz_xyyzz[k];

                g_z_0_zzzzzz_xyzz[k] = -g_zzzzz_xyzz[k] - g_z_0_zzzzz_xyzz[k] * cd_z[k] + g_z_0_zzzzz_xyzzz[k];

                g_z_0_zzzzzz_xzzz[k] = -g_zzzzz_xzzz[k] - g_z_0_zzzzz_xzzz[k] * cd_z[k] + g_z_0_zzzzz_xzzzz[k];

                g_z_0_zzzzzz_yyyy[k] = -g_zzzzz_yyyy[k] - g_z_0_zzzzz_yyyy[k] * cd_z[k] + g_z_0_zzzzz_yyyyz[k];

                g_z_0_zzzzzz_yyyz[k] = -g_zzzzz_yyyz[k] - g_z_0_zzzzz_yyyz[k] * cd_z[k] + g_z_0_zzzzz_yyyzz[k];

                g_z_0_zzzzzz_yyzz[k] = -g_zzzzz_yyzz[k] - g_z_0_zzzzz_yyzz[k] * cd_z[k] + g_z_0_zzzzz_yyzzz[k];

                g_z_0_zzzzzz_yzzz[k] = -g_zzzzz_yzzz[k] - g_z_0_zzzzz_yzzz[k] * cd_z[k] + g_z_0_zzzzz_yzzzz[k];

                g_z_0_zzzzzz_zzzz[k] = -g_zzzzz_zzzz[k] - g_z_0_zzzzz_zzzz[k] * cd_z[k] + g_z_0_zzzzz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

