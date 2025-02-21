#include "ThreeCenterElectronRepulsionGeom010ContrRecXDF.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xdf(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xdf,
                                        const size_t idx_xpf,
                                        const size_t idx_geom_10_xpf,
                                        const size_t idx_geom_10_xpg,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SPF

        const auto pf_off = idx_xpf + i * 30;

        auto g_x_xxx = cbuffer.data(pf_off + 0);

        auto g_x_xxy = cbuffer.data(pf_off + 1);

        auto g_x_xxz = cbuffer.data(pf_off + 2);

        auto g_x_xyy = cbuffer.data(pf_off + 3);

        auto g_x_xyz = cbuffer.data(pf_off + 4);

        auto g_x_xzz = cbuffer.data(pf_off + 5);

        auto g_x_yyy = cbuffer.data(pf_off + 6);

        auto g_x_yyz = cbuffer.data(pf_off + 7);

        auto g_x_yzz = cbuffer.data(pf_off + 8);

        auto g_x_zzz = cbuffer.data(pf_off + 9);

        auto g_y_xxx = cbuffer.data(pf_off + 10);

        auto g_y_xxy = cbuffer.data(pf_off + 11);

        auto g_y_xxz = cbuffer.data(pf_off + 12);

        auto g_y_xyy = cbuffer.data(pf_off + 13);

        auto g_y_xyz = cbuffer.data(pf_off + 14);

        auto g_y_xzz = cbuffer.data(pf_off + 15);

        auto g_y_yyy = cbuffer.data(pf_off + 16);

        auto g_y_yyz = cbuffer.data(pf_off + 17);

        auto g_y_yzz = cbuffer.data(pf_off + 18);

        auto g_y_zzz = cbuffer.data(pf_off + 19);

        auto g_z_xxx = cbuffer.data(pf_off + 20);

        auto g_z_xxy = cbuffer.data(pf_off + 21);

        auto g_z_xxz = cbuffer.data(pf_off + 22);

        auto g_z_xyy = cbuffer.data(pf_off + 23);

        auto g_z_xyz = cbuffer.data(pf_off + 24);

        auto g_z_xzz = cbuffer.data(pf_off + 25);

        auto g_z_yyy = cbuffer.data(pf_off + 26);

        auto g_z_yyz = cbuffer.data(pf_off + 27);

        auto g_z_yzz = cbuffer.data(pf_off + 28);

        auto g_z_zzz = cbuffer.data(pf_off + 29);

        /// Set up components of auxilary buffer : SPF

        const auto pf_geom_10_off = idx_geom_10_xpf + i * 30;

        auto g_x_0_x_xxx = cbuffer.data(pf_geom_10_off + 0 * acomps + 0);

        auto g_x_0_x_xxy = cbuffer.data(pf_geom_10_off + 0 * acomps + 1);

        auto g_x_0_x_xxz = cbuffer.data(pf_geom_10_off + 0 * acomps + 2);

        auto g_x_0_x_xyy = cbuffer.data(pf_geom_10_off + 0 * acomps + 3);

        auto g_x_0_x_xyz = cbuffer.data(pf_geom_10_off + 0 * acomps + 4);

        auto g_x_0_x_xzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 5);

        auto g_x_0_x_yyy = cbuffer.data(pf_geom_10_off + 0 * acomps + 6);

        auto g_x_0_x_yyz = cbuffer.data(pf_geom_10_off + 0 * acomps + 7);

        auto g_x_0_x_yzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 8);

        auto g_x_0_x_zzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 9);

        auto g_x_0_y_xxx = cbuffer.data(pf_geom_10_off + 0 * acomps + 10);

        auto g_x_0_y_xxy = cbuffer.data(pf_geom_10_off + 0 * acomps + 11);

        auto g_x_0_y_xxz = cbuffer.data(pf_geom_10_off + 0 * acomps + 12);

        auto g_x_0_y_xyy = cbuffer.data(pf_geom_10_off + 0 * acomps + 13);

        auto g_x_0_y_xyz = cbuffer.data(pf_geom_10_off + 0 * acomps + 14);

        auto g_x_0_y_xzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 15);

        auto g_x_0_y_yyy = cbuffer.data(pf_geom_10_off + 0 * acomps + 16);

        auto g_x_0_y_yyz = cbuffer.data(pf_geom_10_off + 0 * acomps + 17);

        auto g_x_0_y_yzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 18);

        auto g_x_0_y_zzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 19);

        auto g_x_0_z_xxx = cbuffer.data(pf_geom_10_off + 0 * acomps + 20);

        auto g_x_0_z_xxy = cbuffer.data(pf_geom_10_off + 0 * acomps + 21);

        auto g_x_0_z_xxz = cbuffer.data(pf_geom_10_off + 0 * acomps + 22);

        auto g_x_0_z_xyy = cbuffer.data(pf_geom_10_off + 0 * acomps + 23);

        auto g_x_0_z_xyz = cbuffer.data(pf_geom_10_off + 0 * acomps + 24);

        auto g_x_0_z_xzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 25);

        auto g_x_0_z_yyy = cbuffer.data(pf_geom_10_off + 0 * acomps + 26);

        auto g_x_0_z_yyz = cbuffer.data(pf_geom_10_off + 0 * acomps + 27);

        auto g_x_0_z_yzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 28);

        auto g_x_0_z_zzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 29);

        auto g_y_0_x_xxx = cbuffer.data(pf_geom_10_off + 30 * acomps + 0);

        auto g_y_0_x_xxy = cbuffer.data(pf_geom_10_off + 30 * acomps + 1);

        auto g_y_0_x_xxz = cbuffer.data(pf_geom_10_off + 30 * acomps + 2);

        auto g_y_0_x_xyy = cbuffer.data(pf_geom_10_off + 30 * acomps + 3);

        auto g_y_0_x_xyz = cbuffer.data(pf_geom_10_off + 30 * acomps + 4);

        auto g_y_0_x_xzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 5);

        auto g_y_0_x_yyy = cbuffer.data(pf_geom_10_off + 30 * acomps + 6);

        auto g_y_0_x_yyz = cbuffer.data(pf_geom_10_off + 30 * acomps + 7);

        auto g_y_0_x_yzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 8);

        auto g_y_0_x_zzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 9);

        auto g_y_0_y_xxx = cbuffer.data(pf_geom_10_off + 30 * acomps + 10);

        auto g_y_0_y_xxy = cbuffer.data(pf_geom_10_off + 30 * acomps + 11);

        auto g_y_0_y_xxz = cbuffer.data(pf_geom_10_off + 30 * acomps + 12);

        auto g_y_0_y_xyy = cbuffer.data(pf_geom_10_off + 30 * acomps + 13);

        auto g_y_0_y_xyz = cbuffer.data(pf_geom_10_off + 30 * acomps + 14);

        auto g_y_0_y_xzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 15);

        auto g_y_0_y_yyy = cbuffer.data(pf_geom_10_off + 30 * acomps + 16);

        auto g_y_0_y_yyz = cbuffer.data(pf_geom_10_off + 30 * acomps + 17);

        auto g_y_0_y_yzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 18);

        auto g_y_0_y_zzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 19);

        auto g_y_0_z_xxx = cbuffer.data(pf_geom_10_off + 30 * acomps + 20);

        auto g_y_0_z_xxy = cbuffer.data(pf_geom_10_off + 30 * acomps + 21);

        auto g_y_0_z_xxz = cbuffer.data(pf_geom_10_off + 30 * acomps + 22);

        auto g_y_0_z_xyy = cbuffer.data(pf_geom_10_off + 30 * acomps + 23);

        auto g_y_0_z_xyz = cbuffer.data(pf_geom_10_off + 30 * acomps + 24);

        auto g_y_0_z_xzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 25);

        auto g_y_0_z_yyy = cbuffer.data(pf_geom_10_off + 30 * acomps + 26);

        auto g_y_0_z_yyz = cbuffer.data(pf_geom_10_off + 30 * acomps + 27);

        auto g_y_0_z_yzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 28);

        auto g_y_0_z_zzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 29);

        auto g_z_0_x_xxx = cbuffer.data(pf_geom_10_off + 60 * acomps + 0);

        auto g_z_0_x_xxy = cbuffer.data(pf_geom_10_off + 60 * acomps + 1);

        auto g_z_0_x_xxz = cbuffer.data(pf_geom_10_off + 60 * acomps + 2);

        auto g_z_0_x_xyy = cbuffer.data(pf_geom_10_off + 60 * acomps + 3);

        auto g_z_0_x_xyz = cbuffer.data(pf_geom_10_off + 60 * acomps + 4);

        auto g_z_0_x_xzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 5);

        auto g_z_0_x_yyy = cbuffer.data(pf_geom_10_off + 60 * acomps + 6);

        auto g_z_0_x_yyz = cbuffer.data(pf_geom_10_off + 60 * acomps + 7);

        auto g_z_0_x_yzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 8);

        auto g_z_0_x_zzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 9);

        auto g_z_0_y_xxx = cbuffer.data(pf_geom_10_off + 60 * acomps + 10);

        auto g_z_0_y_xxy = cbuffer.data(pf_geom_10_off + 60 * acomps + 11);

        auto g_z_0_y_xxz = cbuffer.data(pf_geom_10_off + 60 * acomps + 12);

        auto g_z_0_y_xyy = cbuffer.data(pf_geom_10_off + 60 * acomps + 13);

        auto g_z_0_y_xyz = cbuffer.data(pf_geom_10_off + 60 * acomps + 14);

        auto g_z_0_y_xzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 15);

        auto g_z_0_y_yyy = cbuffer.data(pf_geom_10_off + 60 * acomps + 16);

        auto g_z_0_y_yyz = cbuffer.data(pf_geom_10_off + 60 * acomps + 17);

        auto g_z_0_y_yzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 18);

        auto g_z_0_y_zzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 19);

        auto g_z_0_z_xxx = cbuffer.data(pf_geom_10_off + 60 * acomps + 20);

        auto g_z_0_z_xxy = cbuffer.data(pf_geom_10_off + 60 * acomps + 21);

        auto g_z_0_z_xxz = cbuffer.data(pf_geom_10_off + 60 * acomps + 22);

        auto g_z_0_z_xyy = cbuffer.data(pf_geom_10_off + 60 * acomps + 23);

        auto g_z_0_z_xyz = cbuffer.data(pf_geom_10_off + 60 * acomps + 24);

        auto g_z_0_z_xzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 25);

        auto g_z_0_z_yyy = cbuffer.data(pf_geom_10_off + 60 * acomps + 26);

        auto g_z_0_z_yyz = cbuffer.data(pf_geom_10_off + 60 * acomps + 27);

        auto g_z_0_z_yzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 28);

        auto g_z_0_z_zzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 29);

        /// Set up components of auxilary buffer : SPG

        const auto pg_geom_10_off = idx_geom_10_xpg + i * 45;

        auto g_x_0_x_xxxx = cbuffer.data(pg_geom_10_off + 0 * acomps + 0);

        auto g_x_0_x_xxxy = cbuffer.data(pg_geom_10_off + 0 * acomps + 1);

        auto g_x_0_x_xxxz = cbuffer.data(pg_geom_10_off + 0 * acomps + 2);

        auto g_x_0_x_xxyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 3);

        auto g_x_0_x_xxyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 4);

        auto g_x_0_x_xxzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 5);

        auto g_x_0_x_xyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 6);

        auto g_x_0_x_xyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 7);

        auto g_x_0_x_xyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 8);

        auto g_x_0_x_xzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 9);

        auto g_x_0_x_yyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 10);

        auto g_x_0_x_yyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 11);

        auto g_x_0_x_yyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 12);

        auto g_x_0_x_yzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 13);

        auto g_x_0_x_zzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 14);

        auto g_x_0_y_xxxx = cbuffer.data(pg_geom_10_off + 0 * acomps + 15);

        auto g_x_0_y_xxxy = cbuffer.data(pg_geom_10_off + 0 * acomps + 16);

        auto g_x_0_y_xxxz = cbuffer.data(pg_geom_10_off + 0 * acomps + 17);

        auto g_x_0_y_xxyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 18);

        auto g_x_0_y_xxyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 19);

        auto g_x_0_y_xxzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 20);

        auto g_x_0_y_xyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 21);

        auto g_x_0_y_xyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 22);

        auto g_x_0_y_xyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 23);

        auto g_x_0_y_xzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 24);

        auto g_x_0_y_yyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 25);

        auto g_x_0_y_yyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 26);

        auto g_x_0_y_yyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 27);

        auto g_x_0_y_yzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 28);

        auto g_x_0_y_zzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 29);

        auto g_x_0_z_xxxx = cbuffer.data(pg_geom_10_off + 0 * acomps + 30);

        auto g_x_0_z_xxxy = cbuffer.data(pg_geom_10_off + 0 * acomps + 31);

        auto g_x_0_z_xxxz = cbuffer.data(pg_geom_10_off + 0 * acomps + 32);

        auto g_x_0_z_xxyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 33);

        auto g_x_0_z_xxyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 34);

        auto g_x_0_z_xxzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 35);

        auto g_x_0_z_xyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 36);

        auto g_x_0_z_xyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 37);

        auto g_x_0_z_xyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 38);

        auto g_x_0_z_xzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 39);

        auto g_x_0_z_yyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 40);

        auto g_x_0_z_yyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 41);

        auto g_x_0_z_yyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 42);

        auto g_x_0_z_yzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 43);

        auto g_x_0_z_zzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 44);

        auto g_y_0_x_xxxx = cbuffer.data(pg_geom_10_off + 45 * acomps + 0);

        auto g_y_0_x_xxxy = cbuffer.data(pg_geom_10_off + 45 * acomps + 1);

        auto g_y_0_x_xxxz = cbuffer.data(pg_geom_10_off + 45 * acomps + 2);

        auto g_y_0_x_xxyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 3);

        auto g_y_0_x_xxyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 4);

        auto g_y_0_x_xxzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 5);

        auto g_y_0_x_xyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 6);

        auto g_y_0_x_xyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 7);

        auto g_y_0_x_xyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 8);

        auto g_y_0_x_xzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 9);

        auto g_y_0_x_yyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 10);

        auto g_y_0_x_yyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 11);

        auto g_y_0_x_yyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 12);

        auto g_y_0_x_yzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 13);

        auto g_y_0_x_zzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 14);

        auto g_y_0_y_xxxx = cbuffer.data(pg_geom_10_off + 45 * acomps + 15);

        auto g_y_0_y_xxxy = cbuffer.data(pg_geom_10_off + 45 * acomps + 16);

        auto g_y_0_y_xxxz = cbuffer.data(pg_geom_10_off + 45 * acomps + 17);

        auto g_y_0_y_xxyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 18);

        auto g_y_0_y_xxyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 19);

        auto g_y_0_y_xxzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 20);

        auto g_y_0_y_xyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 21);

        auto g_y_0_y_xyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 22);

        auto g_y_0_y_xyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 23);

        auto g_y_0_y_xzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 24);

        auto g_y_0_y_yyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 25);

        auto g_y_0_y_yyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 26);

        auto g_y_0_y_yyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 27);

        auto g_y_0_y_yzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 28);

        auto g_y_0_y_zzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 29);

        auto g_y_0_z_xxxx = cbuffer.data(pg_geom_10_off + 45 * acomps + 30);

        auto g_y_0_z_xxxy = cbuffer.data(pg_geom_10_off + 45 * acomps + 31);

        auto g_y_0_z_xxxz = cbuffer.data(pg_geom_10_off + 45 * acomps + 32);

        auto g_y_0_z_xxyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 33);

        auto g_y_0_z_xxyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 34);

        auto g_y_0_z_xxzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 35);

        auto g_y_0_z_xyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 36);

        auto g_y_0_z_xyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 37);

        auto g_y_0_z_xyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 38);

        auto g_y_0_z_xzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 39);

        auto g_y_0_z_yyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 40);

        auto g_y_0_z_yyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 41);

        auto g_y_0_z_yyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 42);

        auto g_y_0_z_yzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 43);

        auto g_y_0_z_zzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 44);

        auto g_z_0_x_xxxx = cbuffer.data(pg_geom_10_off + 90 * acomps + 0);

        auto g_z_0_x_xxxy = cbuffer.data(pg_geom_10_off + 90 * acomps + 1);

        auto g_z_0_x_xxxz = cbuffer.data(pg_geom_10_off + 90 * acomps + 2);

        auto g_z_0_x_xxyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 3);

        auto g_z_0_x_xxyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 4);

        auto g_z_0_x_xxzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 5);

        auto g_z_0_x_xyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 6);

        auto g_z_0_x_xyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 7);

        auto g_z_0_x_xyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 8);

        auto g_z_0_x_xzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 9);

        auto g_z_0_x_yyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 10);

        auto g_z_0_x_yyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 11);

        auto g_z_0_x_yyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 12);

        auto g_z_0_x_yzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 13);

        auto g_z_0_x_zzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 14);

        auto g_z_0_y_xxxx = cbuffer.data(pg_geom_10_off + 90 * acomps + 15);

        auto g_z_0_y_xxxy = cbuffer.data(pg_geom_10_off + 90 * acomps + 16);

        auto g_z_0_y_xxxz = cbuffer.data(pg_geom_10_off + 90 * acomps + 17);

        auto g_z_0_y_xxyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 18);

        auto g_z_0_y_xxyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 19);

        auto g_z_0_y_xxzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 20);

        auto g_z_0_y_xyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 21);

        auto g_z_0_y_xyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 22);

        auto g_z_0_y_xyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 23);

        auto g_z_0_y_xzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 24);

        auto g_z_0_y_yyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 25);

        auto g_z_0_y_yyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 26);

        auto g_z_0_y_yyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 27);

        auto g_z_0_y_yzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 28);

        auto g_z_0_y_zzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 29);

        auto g_z_0_z_xxxx = cbuffer.data(pg_geom_10_off + 90 * acomps + 30);

        auto g_z_0_z_xxxy = cbuffer.data(pg_geom_10_off + 90 * acomps + 31);

        auto g_z_0_z_xxxz = cbuffer.data(pg_geom_10_off + 90 * acomps + 32);

        auto g_z_0_z_xxyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 33);

        auto g_z_0_z_xxyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 34);

        auto g_z_0_z_xxzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 35);

        auto g_z_0_z_xyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 36);

        auto g_z_0_z_xyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 37);

        auto g_z_0_z_xyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 38);

        auto g_z_0_z_xzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 39);

        auto g_z_0_z_yyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 40);

        auto g_z_0_z_yyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 41);

        auto g_z_0_z_yyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 42);

        auto g_z_0_z_yzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 43);

        auto g_z_0_z_zzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 44);

        /// set up bra offset for contr_buffer_xxdf

        const auto df_geom_10_off = idx_geom_10_xdf + i * 60;

        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_x_0_xx_xxx = cbuffer.data(df_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xx_xxy = cbuffer.data(df_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xx_xxz = cbuffer.data(df_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_xx_xyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xx_xyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xx_xzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_xx_yyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xx_yyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xx_yzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_xx_zzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_x_0_x_xxx, g_x_0_x_xxxx, g_x_0_x_xxxy, g_x_0_x_xxxz, g_x_0_x_xxy, g_x_0_x_xxyy, g_x_0_x_xxyz, g_x_0_x_xxz, g_x_0_x_xxzz, g_x_0_x_xyy, g_x_0_x_xyyy, g_x_0_x_xyyz, g_x_0_x_xyz, g_x_0_x_xyzz, g_x_0_x_xzz, g_x_0_x_xzzz, g_x_0_x_yyy, g_x_0_x_yyz, g_x_0_x_yzz, g_x_0_x_zzz, g_x_0_xx_xxx, g_x_0_xx_xxy, g_x_0_xx_xxz, g_x_0_xx_xyy, g_x_0_xx_xyz, g_x_0_xx_xzz, g_x_0_xx_yyy, g_x_0_xx_yyz, g_x_0_xx_yzz, g_x_0_xx_zzz, g_x_xxx, g_x_xxy, g_x_xxz, g_x_xyy, g_x_xyz, g_x_xzz, g_x_yyy, g_x_yyz, g_x_yzz, g_x_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xx_xxx[k] = -g_x_xxx[k] - g_x_0_x_xxx[k] * cd_x[k] + g_x_0_x_xxxx[k];

            g_x_0_xx_xxy[k] = -g_x_xxy[k] - g_x_0_x_xxy[k] * cd_x[k] + g_x_0_x_xxxy[k];

            g_x_0_xx_xxz[k] = -g_x_xxz[k] - g_x_0_x_xxz[k] * cd_x[k] + g_x_0_x_xxxz[k];

            g_x_0_xx_xyy[k] = -g_x_xyy[k] - g_x_0_x_xyy[k] * cd_x[k] + g_x_0_x_xxyy[k];

            g_x_0_xx_xyz[k] = -g_x_xyz[k] - g_x_0_x_xyz[k] * cd_x[k] + g_x_0_x_xxyz[k];

            g_x_0_xx_xzz[k] = -g_x_xzz[k] - g_x_0_x_xzz[k] * cd_x[k] + g_x_0_x_xxzz[k];

            g_x_0_xx_yyy[k] = -g_x_yyy[k] - g_x_0_x_yyy[k] * cd_x[k] + g_x_0_x_xyyy[k];

            g_x_0_xx_yyz[k] = -g_x_yyz[k] - g_x_0_x_yyz[k] * cd_x[k] + g_x_0_x_xyyz[k];

            g_x_0_xx_yzz[k] = -g_x_yzz[k] - g_x_0_x_yzz[k] * cd_x[k] + g_x_0_x_xyzz[k];

            g_x_0_xx_zzz[k] = -g_x_zzz[k] - g_x_0_x_zzz[k] * cd_x[k] + g_x_0_x_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_x_0_xy_xxx = cbuffer.data(df_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xy_xxy = cbuffer.data(df_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_xy_xxz = cbuffer.data(df_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xy_xyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xy_xyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_xy_xzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xy_yyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xy_yyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_xy_yzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_xy_zzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 19);

        #pragma omp simd aligned(cd_y, g_x_0_x_xxx, g_x_0_x_xxxy, g_x_0_x_xxy, g_x_0_x_xxyy, g_x_0_x_xxyz, g_x_0_x_xxz, g_x_0_x_xyy, g_x_0_x_xyyy, g_x_0_x_xyyz, g_x_0_x_xyz, g_x_0_x_xyzz, g_x_0_x_xzz, g_x_0_x_yyy, g_x_0_x_yyyy, g_x_0_x_yyyz, g_x_0_x_yyz, g_x_0_x_yyzz, g_x_0_x_yzz, g_x_0_x_yzzz, g_x_0_x_zzz, g_x_0_xy_xxx, g_x_0_xy_xxy, g_x_0_xy_xxz, g_x_0_xy_xyy, g_x_0_xy_xyz, g_x_0_xy_xzz, g_x_0_xy_yyy, g_x_0_xy_yyz, g_x_0_xy_yzz, g_x_0_xy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xy_xxx[k] = -g_x_0_x_xxx[k] * cd_y[k] + g_x_0_x_xxxy[k];

            g_x_0_xy_xxy[k] = -g_x_0_x_xxy[k] * cd_y[k] + g_x_0_x_xxyy[k];

            g_x_0_xy_xxz[k] = -g_x_0_x_xxz[k] * cd_y[k] + g_x_0_x_xxyz[k];

            g_x_0_xy_xyy[k] = -g_x_0_x_xyy[k] * cd_y[k] + g_x_0_x_xyyy[k];

            g_x_0_xy_xyz[k] = -g_x_0_x_xyz[k] * cd_y[k] + g_x_0_x_xyyz[k];

            g_x_0_xy_xzz[k] = -g_x_0_x_xzz[k] * cd_y[k] + g_x_0_x_xyzz[k];

            g_x_0_xy_yyy[k] = -g_x_0_x_yyy[k] * cd_y[k] + g_x_0_x_yyyy[k];

            g_x_0_xy_yyz[k] = -g_x_0_x_yyz[k] * cd_y[k] + g_x_0_x_yyyz[k];

            g_x_0_xy_yzz[k] = -g_x_0_x_yzz[k] * cd_y[k] + g_x_0_x_yyzz[k];

            g_x_0_xy_zzz[k] = -g_x_0_x_zzz[k] * cd_y[k] + g_x_0_x_yzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_x_0_xz_xxx = cbuffer.data(df_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_xz_xxy = cbuffer.data(df_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_xz_xxz = cbuffer.data(df_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_xz_xyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_xz_xyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_xz_xzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_xz_yyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_xz_yyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_xz_yzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_xz_zzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 29);

        #pragma omp simd aligned(cd_z, g_x_0_x_xxx, g_x_0_x_xxxz, g_x_0_x_xxy, g_x_0_x_xxyz, g_x_0_x_xxz, g_x_0_x_xxzz, g_x_0_x_xyy, g_x_0_x_xyyz, g_x_0_x_xyz, g_x_0_x_xyzz, g_x_0_x_xzz, g_x_0_x_xzzz, g_x_0_x_yyy, g_x_0_x_yyyz, g_x_0_x_yyz, g_x_0_x_yyzz, g_x_0_x_yzz, g_x_0_x_yzzz, g_x_0_x_zzz, g_x_0_x_zzzz, g_x_0_xz_xxx, g_x_0_xz_xxy, g_x_0_xz_xxz, g_x_0_xz_xyy, g_x_0_xz_xyz, g_x_0_xz_xzz, g_x_0_xz_yyy, g_x_0_xz_yyz, g_x_0_xz_yzz, g_x_0_xz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xz_xxx[k] = -g_x_0_x_xxx[k] * cd_z[k] + g_x_0_x_xxxz[k];

            g_x_0_xz_xxy[k] = -g_x_0_x_xxy[k] * cd_z[k] + g_x_0_x_xxyz[k];

            g_x_0_xz_xxz[k] = -g_x_0_x_xxz[k] * cd_z[k] + g_x_0_x_xxzz[k];

            g_x_0_xz_xyy[k] = -g_x_0_x_xyy[k] * cd_z[k] + g_x_0_x_xyyz[k];

            g_x_0_xz_xyz[k] = -g_x_0_x_xyz[k] * cd_z[k] + g_x_0_x_xyzz[k];

            g_x_0_xz_xzz[k] = -g_x_0_x_xzz[k] * cd_z[k] + g_x_0_x_xzzz[k];

            g_x_0_xz_yyy[k] = -g_x_0_x_yyy[k] * cd_z[k] + g_x_0_x_yyyz[k];

            g_x_0_xz_yyz[k] = -g_x_0_x_yyz[k] * cd_z[k] + g_x_0_x_yyzz[k];

            g_x_0_xz_yzz[k] = -g_x_0_x_yzz[k] * cd_z[k] + g_x_0_x_yzzz[k];

            g_x_0_xz_zzz[k] = -g_x_0_x_zzz[k] * cd_z[k] + g_x_0_x_zzzz[k];
        }

        /// Set up 30-40 components of targeted buffer : cbuffer.data(

        auto g_x_0_yy_xxx = cbuffer.data(df_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_yy_xxy = cbuffer.data(df_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_yy_xxz = cbuffer.data(df_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_yy_xyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_yy_xyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_yy_xzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 35);

        auto g_x_0_yy_yyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_yy_yyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_yy_yzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_yy_zzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 39);

        #pragma omp simd aligned(cd_y, g_x_0_y_xxx, g_x_0_y_xxxy, g_x_0_y_xxy, g_x_0_y_xxyy, g_x_0_y_xxyz, g_x_0_y_xxz, g_x_0_y_xyy, g_x_0_y_xyyy, g_x_0_y_xyyz, g_x_0_y_xyz, g_x_0_y_xyzz, g_x_0_y_xzz, g_x_0_y_yyy, g_x_0_y_yyyy, g_x_0_y_yyyz, g_x_0_y_yyz, g_x_0_y_yyzz, g_x_0_y_yzz, g_x_0_y_yzzz, g_x_0_y_zzz, g_x_0_yy_xxx, g_x_0_yy_xxy, g_x_0_yy_xxz, g_x_0_yy_xyy, g_x_0_yy_xyz, g_x_0_yy_xzz, g_x_0_yy_yyy, g_x_0_yy_yyz, g_x_0_yy_yzz, g_x_0_yy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yy_xxx[k] = -g_x_0_y_xxx[k] * cd_y[k] + g_x_0_y_xxxy[k];

            g_x_0_yy_xxy[k] = -g_x_0_y_xxy[k] * cd_y[k] + g_x_0_y_xxyy[k];

            g_x_0_yy_xxz[k] = -g_x_0_y_xxz[k] * cd_y[k] + g_x_0_y_xxyz[k];

            g_x_0_yy_xyy[k] = -g_x_0_y_xyy[k] * cd_y[k] + g_x_0_y_xyyy[k];

            g_x_0_yy_xyz[k] = -g_x_0_y_xyz[k] * cd_y[k] + g_x_0_y_xyyz[k];

            g_x_0_yy_xzz[k] = -g_x_0_y_xzz[k] * cd_y[k] + g_x_0_y_xyzz[k];

            g_x_0_yy_yyy[k] = -g_x_0_y_yyy[k] * cd_y[k] + g_x_0_y_yyyy[k];

            g_x_0_yy_yyz[k] = -g_x_0_y_yyz[k] * cd_y[k] + g_x_0_y_yyyz[k];

            g_x_0_yy_yzz[k] = -g_x_0_y_yzz[k] * cd_y[k] + g_x_0_y_yyzz[k];

            g_x_0_yy_zzz[k] = -g_x_0_y_zzz[k] * cd_y[k] + g_x_0_y_yzzz[k];
        }

        /// Set up 40-50 components of targeted buffer : cbuffer.data(

        auto g_x_0_yz_xxx = cbuffer.data(df_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_yz_xxy = cbuffer.data(df_geom_10_off + 0 * acomps  + 41);

        auto g_x_0_yz_xxz = cbuffer.data(df_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_yz_xyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_yz_xyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 44);

        auto g_x_0_yz_xzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_yz_yyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_yz_yyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 47);

        auto g_x_0_yz_yzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_yz_zzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 49);

        #pragma omp simd aligned(cd_y, g_x_0_yz_xxx, g_x_0_yz_xxy, g_x_0_yz_xxz, g_x_0_yz_xyy, g_x_0_yz_xyz, g_x_0_yz_xzz, g_x_0_yz_yyy, g_x_0_yz_yyz, g_x_0_yz_yzz, g_x_0_yz_zzz, g_x_0_z_xxx, g_x_0_z_xxxy, g_x_0_z_xxy, g_x_0_z_xxyy, g_x_0_z_xxyz, g_x_0_z_xxz, g_x_0_z_xyy, g_x_0_z_xyyy, g_x_0_z_xyyz, g_x_0_z_xyz, g_x_0_z_xyzz, g_x_0_z_xzz, g_x_0_z_yyy, g_x_0_z_yyyy, g_x_0_z_yyyz, g_x_0_z_yyz, g_x_0_z_yyzz, g_x_0_z_yzz, g_x_0_z_yzzz, g_x_0_z_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yz_xxx[k] = -g_x_0_z_xxx[k] * cd_y[k] + g_x_0_z_xxxy[k];

            g_x_0_yz_xxy[k] = -g_x_0_z_xxy[k] * cd_y[k] + g_x_0_z_xxyy[k];

            g_x_0_yz_xxz[k] = -g_x_0_z_xxz[k] * cd_y[k] + g_x_0_z_xxyz[k];

            g_x_0_yz_xyy[k] = -g_x_0_z_xyy[k] * cd_y[k] + g_x_0_z_xyyy[k];

            g_x_0_yz_xyz[k] = -g_x_0_z_xyz[k] * cd_y[k] + g_x_0_z_xyyz[k];

            g_x_0_yz_xzz[k] = -g_x_0_z_xzz[k] * cd_y[k] + g_x_0_z_xyzz[k];

            g_x_0_yz_yyy[k] = -g_x_0_z_yyy[k] * cd_y[k] + g_x_0_z_yyyy[k];

            g_x_0_yz_yyz[k] = -g_x_0_z_yyz[k] * cd_y[k] + g_x_0_z_yyyz[k];

            g_x_0_yz_yzz[k] = -g_x_0_z_yzz[k] * cd_y[k] + g_x_0_z_yyzz[k];

            g_x_0_yz_zzz[k] = -g_x_0_z_zzz[k] * cd_y[k] + g_x_0_z_yzzz[k];
        }

        /// Set up 50-60 components of targeted buffer : cbuffer.data(

        auto g_x_0_zz_xxx = cbuffer.data(df_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_zz_xxy = cbuffer.data(df_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_zz_xxz = cbuffer.data(df_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_zz_xyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 53);

        auto g_x_0_zz_xyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_zz_xzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_zz_yyy = cbuffer.data(df_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_zz_yyz = cbuffer.data(df_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_zz_yzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_zz_zzz = cbuffer.data(df_geom_10_off + 0 * acomps  + 59);

        #pragma omp simd aligned(cd_z, g_x_0_z_xxx, g_x_0_z_xxxz, g_x_0_z_xxy, g_x_0_z_xxyz, g_x_0_z_xxz, g_x_0_z_xxzz, g_x_0_z_xyy, g_x_0_z_xyyz, g_x_0_z_xyz, g_x_0_z_xyzz, g_x_0_z_xzz, g_x_0_z_xzzz, g_x_0_z_yyy, g_x_0_z_yyyz, g_x_0_z_yyz, g_x_0_z_yyzz, g_x_0_z_yzz, g_x_0_z_yzzz, g_x_0_z_zzz, g_x_0_z_zzzz, g_x_0_zz_xxx, g_x_0_zz_xxy, g_x_0_zz_xxz, g_x_0_zz_xyy, g_x_0_zz_xyz, g_x_0_zz_xzz, g_x_0_zz_yyy, g_x_0_zz_yyz, g_x_0_zz_yzz, g_x_0_zz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zz_xxx[k] = -g_x_0_z_xxx[k] * cd_z[k] + g_x_0_z_xxxz[k];

            g_x_0_zz_xxy[k] = -g_x_0_z_xxy[k] * cd_z[k] + g_x_0_z_xxyz[k];

            g_x_0_zz_xxz[k] = -g_x_0_z_xxz[k] * cd_z[k] + g_x_0_z_xxzz[k];

            g_x_0_zz_xyy[k] = -g_x_0_z_xyy[k] * cd_z[k] + g_x_0_z_xyyz[k];

            g_x_0_zz_xyz[k] = -g_x_0_z_xyz[k] * cd_z[k] + g_x_0_z_xyzz[k];

            g_x_0_zz_xzz[k] = -g_x_0_z_xzz[k] * cd_z[k] + g_x_0_z_xzzz[k];

            g_x_0_zz_yyy[k] = -g_x_0_z_yyy[k] * cd_z[k] + g_x_0_z_yyyz[k];

            g_x_0_zz_yyz[k] = -g_x_0_z_yyz[k] * cd_z[k] + g_x_0_z_yyzz[k];

            g_x_0_zz_yzz[k] = -g_x_0_z_yzz[k] * cd_z[k] + g_x_0_z_yzzz[k];

            g_x_0_zz_zzz[k] = -g_x_0_z_zzz[k] * cd_z[k] + g_x_0_z_zzzz[k];
        }
        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_y_0_xx_xxx = cbuffer.data(df_geom_10_off + 60 * acomps  + 0);

        auto g_y_0_xx_xxy = cbuffer.data(df_geom_10_off + 60 * acomps  + 1);

        auto g_y_0_xx_xxz = cbuffer.data(df_geom_10_off + 60 * acomps  + 2);

        auto g_y_0_xx_xyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 3);

        auto g_y_0_xx_xyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 4);

        auto g_y_0_xx_xzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 5);

        auto g_y_0_xx_yyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 6);

        auto g_y_0_xx_yyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 7);

        auto g_y_0_xx_yzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 8);

        auto g_y_0_xx_zzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_y_0_x_xxx, g_y_0_x_xxxx, g_y_0_x_xxxy, g_y_0_x_xxxz, g_y_0_x_xxy, g_y_0_x_xxyy, g_y_0_x_xxyz, g_y_0_x_xxz, g_y_0_x_xxzz, g_y_0_x_xyy, g_y_0_x_xyyy, g_y_0_x_xyyz, g_y_0_x_xyz, g_y_0_x_xyzz, g_y_0_x_xzz, g_y_0_x_xzzz, g_y_0_x_yyy, g_y_0_x_yyz, g_y_0_x_yzz, g_y_0_x_zzz, g_y_0_xx_xxx, g_y_0_xx_xxy, g_y_0_xx_xxz, g_y_0_xx_xyy, g_y_0_xx_xyz, g_y_0_xx_xzz, g_y_0_xx_yyy, g_y_0_xx_yyz, g_y_0_xx_yzz, g_y_0_xx_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xx_xxx[k] = -g_y_0_x_xxx[k] * cd_x[k] + g_y_0_x_xxxx[k];

            g_y_0_xx_xxy[k] = -g_y_0_x_xxy[k] * cd_x[k] + g_y_0_x_xxxy[k];

            g_y_0_xx_xxz[k] = -g_y_0_x_xxz[k] * cd_x[k] + g_y_0_x_xxxz[k];

            g_y_0_xx_xyy[k] = -g_y_0_x_xyy[k] * cd_x[k] + g_y_0_x_xxyy[k];

            g_y_0_xx_xyz[k] = -g_y_0_x_xyz[k] * cd_x[k] + g_y_0_x_xxyz[k];

            g_y_0_xx_xzz[k] = -g_y_0_x_xzz[k] * cd_x[k] + g_y_0_x_xxzz[k];

            g_y_0_xx_yyy[k] = -g_y_0_x_yyy[k] * cd_x[k] + g_y_0_x_xyyy[k];

            g_y_0_xx_yyz[k] = -g_y_0_x_yyz[k] * cd_x[k] + g_y_0_x_xyyz[k];

            g_y_0_xx_yzz[k] = -g_y_0_x_yzz[k] * cd_x[k] + g_y_0_x_xyzz[k];

            g_y_0_xx_zzz[k] = -g_y_0_x_zzz[k] * cd_x[k] + g_y_0_x_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_y_0_xy_xxx = cbuffer.data(df_geom_10_off + 60 * acomps  + 10);

        auto g_y_0_xy_xxy = cbuffer.data(df_geom_10_off + 60 * acomps  + 11);

        auto g_y_0_xy_xxz = cbuffer.data(df_geom_10_off + 60 * acomps  + 12);

        auto g_y_0_xy_xyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 13);

        auto g_y_0_xy_xyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 14);

        auto g_y_0_xy_xzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 15);

        auto g_y_0_xy_yyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 16);

        auto g_y_0_xy_yyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 17);

        auto g_y_0_xy_yzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 18);

        auto g_y_0_xy_zzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 19);

        #pragma omp simd aligned(cd_x, g_y_0_xy_xxx, g_y_0_xy_xxy, g_y_0_xy_xxz, g_y_0_xy_xyy, g_y_0_xy_xyz, g_y_0_xy_xzz, g_y_0_xy_yyy, g_y_0_xy_yyz, g_y_0_xy_yzz, g_y_0_xy_zzz, g_y_0_y_xxx, g_y_0_y_xxxx, g_y_0_y_xxxy, g_y_0_y_xxxz, g_y_0_y_xxy, g_y_0_y_xxyy, g_y_0_y_xxyz, g_y_0_y_xxz, g_y_0_y_xxzz, g_y_0_y_xyy, g_y_0_y_xyyy, g_y_0_y_xyyz, g_y_0_y_xyz, g_y_0_y_xyzz, g_y_0_y_xzz, g_y_0_y_xzzz, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yzz, g_y_0_y_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xy_xxx[k] = -g_y_0_y_xxx[k] * cd_x[k] + g_y_0_y_xxxx[k];

            g_y_0_xy_xxy[k] = -g_y_0_y_xxy[k] * cd_x[k] + g_y_0_y_xxxy[k];

            g_y_0_xy_xxz[k] = -g_y_0_y_xxz[k] * cd_x[k] + g_y_0_y_xxxz[k];

            g_y_0_xy_xyy[k] = -g_y_0_y_xyy[k] * cd_x[k] + g_y_0_y_xxyy[k];

            g_y_0_xy_xyz[k] = -g_y_0_y_xyz[k] * cd_x[k] + g_y_0_y_xxyz[k];

            g_y_0_xy_xzz[k] = -g_y_0_y_xzz[k] * cd_x[k] + g_y_0_y_xxzz[k];

            g_y_0_xy_yyy[k] = -g_y_0_y_yyy[k] * cd_x[k] + g_y_0_y_xyyy[k];

            g_y_0_xy_yyz[k] = -g_y_0_y_yyz[k] * cd_x[k] + g_y_0_y_xyyz[k];

            g_y_0_xy_yzz[k] = -g_y_0_y_yzz[k] * cd_x[k] + g_y_0_y_xyzz[k];

            g_y_0_xy_zzz[k] = -g_y_0_y_zzz[k] * cd_x[k] + g_y_0_y_xzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_y_0_xz_xxx = cbuffer.data(df_geom_10_off + 60 * acomps  + 20);

        auto g_y_0_xz_xxy = cbuffer.data(df_geom_10_off + 60 * acomps  + 21);

        auto g_y_0_xz_xxz = cbuffer.data(df_geom_10_off + 60 * acomps  + 22);

        auto g_y_0_xz_xyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 23);

        auto g_y_0_xz_xyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 24);

        auto g_y_0_xz_xzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 25);

        auto g_y_0_xz_yyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 26);

        auto g_y_0_xz_yyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 27);

        auto g_y_0_xz_yzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 28);

        auto g_y_0_xz_zzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_y_0_xz_xxx, g_y_0_xz_xxy, g_y_0_xz_xxz, g_y_0_xz_xyy, g_y_0_xz_xyz, g_y_0_xz_xzz, g_y_0_xz_yyy, g_y_0_xz_yyz, g_y_0_xz_yzz, g_y_0_xz_zzz, g_y_0_z_xxx, g_y_0_z_xxxx, g_y_0_z_xxxy, g_y_0_z_xxxz, g_y_0_z_xxy, g_y_0_z_xxyy, g_y_0_z_xxyz, g_y_0_z_xxz, g_y_0_z_xxzz, g_y_0_z_xyy, g_y_0_z_xyyy, g_y_0_z_xyyz, g_y_0_z_xyz, g_y_0_z_xyzz, g_y_0_z_xzz, g_y_0_z_xzzz, g_y_0_z_yyy, g_y_0_z_yyz, g_y_0_z_yzz, g_y_0_z_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xz_xxx[k] = -g_y_0_z_xxx[k] * cd_x[k] + g_y_0_z_xxxx[k];

            g_y_0_xz_xxy[k] = -g_y_0_z_xxy[k] * cd_x[k] + g_y_0_z_xxxy[k];

            g_y_0_xz_xxz[k] = -g_y_0_z_xxz[k] * cd_x[k] + g_y_0_z_xxxz[k];

            g_y_0_xz_xyy[k] = -g_y_0_z_xyy[k] * cd_x[k] + g_y_0_z_xxyy[k];

            g_y_0_xz_xyz[k] = -g_y_0_z_xyz[k] * cd_x[k] + g_y_0_z_xxyz[k];

            g_y_0_xz_xzz[k] = -g_y_0_z_xzz[k] * cd_x[k] + g_y_0_z_xxzz[k];

            g_y_0_xz_yyy[k] = -g_y_0_z_yyy[k] * cd_x[k] + g_y_0_z_xyyy[k];

            g_y_0_xz_yyz[k] = -g_y_0_z_yyz[k] * cd_x[k] + g_y_0_z_xyyz[k];

            g_y_0_xz_yzz[k] = -g_y_0_z_yzz[k] * cd_x[k] + g_y_0_z_xyzz[k];

            g_y_0_xz_zzz[k] = -g_y_0_z_zzz[k] * cd_x[k] + g_y_0_z_xzzz[k];
        }

        /// Set up 30-40 components of targeted buffer : cbuffer.data(

        auto g_y_0_yy_xxx = cbuffer.data(df_geom_10_off + 60 * acomps  + 30);

        auto g_y_0_yy_xxy = cbuffer.data(df_geom_10_off + 60 * acomps  + 31);

        auto g_y_0_yy_xxz = cbuffer.data(df_geom_10_off + 60 * acomps  + 32);

        auto g_y_0_yy_xyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 33);

        auto g_y_0_yy_xyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 34);

        auto g_y_0_yy_xzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 35);

        auto g_y_0_yy_yyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 36);

        auto g_y_0_yy_yyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 37);

        auto g_y_0_yy_yzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 38);

        auto g_y_0_yy_zzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 39);

        #pragma omp simd aligned(cd_y, g_y_0_y_xxx, g_y_0_y_xxxy, g_y_0_y_xxy, g_y_0_y_xxyy, g_y_0_y_xxyz, g_y_0_y_xxz, g_y_0_y_xyy, g_y_0_y_xyyy, g_y_0_y_xyyz, g_y_0_y_xyz, g_y_0_y_xyzz, g_y_0_y_xzz, g_y_0_y_yyy, g_y_0_y_yyyy, g_y_0_y_yyyz, g_y_0_y_yyz, g_y_0_y_yyzz, g_y_0_y_yzz, g_y_0_y_yzzz, g_y_0_y_zzz, g_y_0_yy_xxx, g_y_0_yy_xxy, g_y_0_yy_xxz, g_y_0_yy_xyy, g_y_0_yy_xyz, g_y_0_yy_xzz, g_y_0_yy_yyy, g_y_0_yy_yyz, g_y_0_yy_yzz, g_y_0_yy_zzz, g_y_xxx, g_y_xxy, g_y_xxz, g_y_xyy, g_y_xyz, g_y_xzz, g_y_yyy, g_y_yyz, g_y_yzz, g_y_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yy_xxx[k] = -g_y_xxx[k] - g_y_0_y_xxx[k] * cd_y[k] + g_y_0_y_xxxy[k];

            g_y_0_yy_xxy[k] = -g_y_xxy[k] - g_y_0_y_xxy[k] * cd_y[k] + g_y_0_y_xxyy[k];

            g_y_0_yy_xxz[k] = -g_y_xxz[k] - g_y_0_y_xxz[k] * cd_y[k] + g_y_0_y_xxyz[k];

            g_y_0_yy_xyy[k] = -g_y_xyy[k] - g_y_0_y_xyy[k] * cd_y[k] + g_y_0_y_xyyy[k];

            g_y_0_yy_xyz[k] = -g_y_xyz[k] - g_y_0_y_xyz[k] * cd_y[k] + g_y_0_y_xyyz[k];

            g_y_0_yy_xzz[k] = -g_y_xzz[k] - g_y_0_y_xzz[k] * cd_y[k] + g_y_0_y_xyzz[k];

            g_y_0_yy_yyy[k] = -g_y_yyy[k] - g_y_0_y_yyy[k] * cd_y[k] + g_y_0_y_yyyy[k];

            g_y_0_yy_yyz[k] = -g_y_yyz[k] - g_y_0_y_yyz[k] * cd_y[k] + g_y_0_y_yyyz[k];

            g_y_0_yy_yzz[k] = -g_y_yzz[k] - g_y_0_y_yzz[k] * cd_y[k] + g_y_0_y_yyzz[k];

            g_y_0_yy_zzz[k] = -g_y_zzz[k] - g_y_0_y_zzz[k] * cd_y[k] + g_y_0_y_yzzz[k];
        }

        /// Set up 40-50 components of targeted buffer : cbuffer.data(

        auto g_y_0_yz_xxx = cbuffer.data(df_geom_10_off + 60 * acomps  + 40);

        auto g_y_0_yz_xxy = cbuffer.data(df_geom_10_off + 60 * acomps  + 41);

        auto g_y_0_yz_xxz = cbuffer.data(df_geom_10_off + 60 * acomps  + 42);

        auto g_y_0_yz_xyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 43);

        auto g_y_0_yz_xyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 44);

        auto g_y_0_yz_xzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 45);

        auto g_y_0_yz_yyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 46);

        auto g_y_0_yz_yyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 47);

        auto g_y_0_yz_yzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 48);

        auto g_y_0_yz_zzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 49);

        #pragma omp simd aligned(cd_z, g_y_0_y_xxx, g_y_0_y_xxxz, g_y_0_y_xxy, g_y_0_y_xxyz, g_y_0_y_xxz, g_y_0_y_xxzz, g_y_0_y_xyy, g_y_0_y_xyyz, g_y_0_y_xyz, g_y_0_y_xyzz, g_y_0_y_xzz, g_y_0_y_xzzz, g_y_0_y_yyy, g_y_0_y_yyyz, g_y_0_y_yyz, g_y_0_y_yyzz, g_y_0_y_yzz, g_y_0_y_yzzz, g_y_0_y_zzz, g_y_0_y_zzzz, g_y_0_yz_xxx, g_y_0_yz_xxy, g_y_0_yz_xxz, g_y_0_yz_xyy, g_y_0_yz_xyz, g_y_0_yz_xzz, g_y_0_yz_yyy, g_y_0_yz_yyz, g_y_0_yz_yzz, g_y_0_yz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yz_xxx[k] = -g_y_0_y_xxx[k] * cd_z[k] + g_y_0_y_xxxz[k];

            g_y_0_yz_xxy[k] = -g_y_0_y_xxy[k] * cd_z[k] + g_y_0_y_xxyz[k];

            g_y_0_yz_xxz[k] = -g_y_0_y_xxz[k] * cd_z[k] + g_y_0_y_xxzz[k];

            g_y_0_yz_xyy[k] = -g_y_0_y_xyy[k] * cd_z[k] + g_y_0_y_xyyz[k];

            g_y_0_yz_xyz[k] = -g_y_0_y_xyz[k] * cd_z[k] + g_y_0_y_xyzz[k];

            g_y_0_yz_xzz[k] = -g_y_0_y_xzz[k] * cd_z[k] + g_y_0_y_xzzz[k];

            g_y_0_yz_yyy[k] = -g_y_0_y_yyy[k] * cd_z[k] + g_y_0_y_yyyz[k];

            g_y_0_yz_yyz[k] = -g_y_0_y_yyz[k] * cd_z[k] + g_y_0_y_yyzz[k];

            g_y_0_yz_yzz[k] = -g_y_0_y_yzz[k] * cd_z[k] + g_y_0_y_yzzz[k];

            g_y_0_yz_zzz[k] = -g_y_0_y_zzz[k] * cd_z[k] + g_y_0_y_zzzz[k];
        }

        /// Set up 50-60 components of targeted buffer : cbuffer.data(

        auto g_y_0_zz_xxx = cbuffer.data(df_geom_10_off + 60 * acomps  + 50);

        auto g_y_0_zz_xxy = cbuffer.data(df_geom_10_off + 60 * acomps  + 51);

        auto g_y_0_zz_xxz = cbuffer.data(df_geom_10_off + 60 * acomps  + 52);

        auto g_y_0_zz_xyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 53);

        auto g_y_0_zz_xyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 54);

        auto g_y_0_zz_xzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 55);

        auto g_y_0_zz_yyy = cbuffer.data(df_geom_10_off + 60 * acomps  + 56);

        auto g_y_0_zz_yyz = cbuffer.data(df_geom_10_off + 60 * acomps  + 57);

        auto g_y_0_zz_yzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 58);

        auto g_y_0_zz_zzz = cbuffer.data(df_geom_10_off + 60 * acomps  + 59);

        #pragma omp simd aligned(cd_z, g_y_0_z_xxx, g_y_0_z_xxxz, g_y_0_z_xxy, g_y_0_z_xxyz, g_y_0_z_xxz, g_y_0_z_xxzz, g_y_0_z_xyy, g_y_0_z_xyyz, g_y_0_z_xyz, g_y_0_z_xyzz, g_y_0_z_xzz, g_y_0_z_xzzz, g_y_0_z_yyy, g_y_0_z_yyyz, g_y_0_z_yyz, g_y_0_z_yyzz, g_y_0_z_yzz, g_y_0_z_yzzz, g_y_0_z_zzz, g_y_0_z_zzzz, g_y_0_zz_xxx, g_y_0_zz_xxy, g_y_0_zz_xxz, g_y_0_zz_xyy, g_y_0_zz_xyz, g_y_0_zz_xzz, g_y_0_zz_yyy, g_y_0_zz_yyz, g_y_0_zz_yzz, g_y_0_zz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zz_xxx[k] = -g_y_0_z_xxx[k] * cd_z[k] + g_y_0_z_xxxz[k];

            g_y_0_zz_xxy[k] = -g_y_0_z_xxy[k] * cd_z[k] + g_y_0_z_xxyz[k];

            g_y_0_zz_xxz[k] = -g_y_0_z_xxz[k] * cd_z[k] + g_y_0_z_xxzz[k];

            g_y_0_zz_xyy[k] = -g_y_0_z_xyy[k] * cd_z[k] + g_y_0_z_xyyz[k];

            g_y_0_zz_xyz[k] = -g_y_0_z_xyz[k] * cd_z[k] + g_y_0_z_xyzz[k];

            g_y_0_zz_xzz[k] = -g_y_0_z_xzz[k] * cd_z[k] + g_y_0_z_xzzz[k];

            g_y_0_zz_yyy[k] = -g_y_0_z_yyy[k] * cd_z[k] + g_y_0_z_yyyz[k];

            g_y_0_zz_yyz[k] = -g_y_0_z_yyz[k] * cd_z[k] + g_y_0_z_yyzz[k];

            g_y_0_zz_yzz[k] = -g_y_0_z_yzz[k] * cd_z[k] + g_y_0_z_yzzz[k];

            g_y_0_zz_zzz[k] = -g_y_0_z_zzz[k] * cd_z[k] + g_y_0_z_zzzz[k];
        }
        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_z_0_xx_xxx = cbuffer.data(df_geom_10_off + 120 * acomps  + 0);

        auto g_z_0_xx_xxy = cbuffer.data(df_geom_10_off + 120 * acomps  + 1);

        auto g_z_0_xx_xxz = cbuffer.data(df_geom_10_off + 120 * acomps  + 2);

        auto g_z_0_xx_xyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 3);

        auto g_z_0_xx_xyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 4);

        auto g_z_0_xx_xzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 5);

        auto g_z_0_xx_yyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 6);

        auto g_z_0_xx_yyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 7);

        auto g_z_0_xx_yzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 8);

        auto g_z_0_xx_zzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_z_0_x_xxx, g_z_0_x_xxxx, g_z_0_x_xxxy, g_z_0_x_xxxz, g_z_0_x_xxy, g_z_0_x_xxyy, g_z_0_x_xxyz, g_z_0_x_xxz, g_z_0_x_xxzz, g_z_0_x_xyy, g_z_0_x_xyyy, g_z_0_x_xyyz, g_z_0_x_xyz, g_z_0_x_xyzz, g_z_0_x_xzz, g_z_0_x_xzzz, g_z_0_x_yyy, g_z_0_x_yyz, g_z_0_x_yzz, g_z_0_x_zzz, g_z_0_xx_xxx, g_z_0_xx_xxy, g_z_0_xx_xxz, g_z_0_xx_xyy, g_z_0_xx_xyz, g_z_0_xx_xzz, g_z_0_xx_yyy, g_z_0_xx_yyz, g_z_0_xx_yzz, g_z_0_xx_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xx_xxx[k] = -g_z_0_x_xxx[k] * cd_x[k] + g_z_0_x_xxxx[k];

            g_z_0_xx_xxy[k] = -g_z_0_x_xxy[k] * cd_x[k] + g_z_0_x_xxxy[k];

            g_z_0_xx_xxz[k] = -g_z_0_x_xxz[k] * cd_x[k] + g_z_0_x_xxxz[k];

            g_z_0_xx_xyy[k] = -g_z_0_x_xyy[k] * cd_x[k] + g_z_0_x_xxyy[k];

            g_z_0_xx_xyz[k] = -g_z_0_x_xyz[k] * cd_x[k] + g_z_0_x_xxyz[k];

            g_z_0_xx_xzz[k] = -g_z_0_x_xzz[k] * cd_x[k] + g_z_0_x_xxzz[k];

            g_z_0_xx_yyy[k] = -g_z_0_x_yyy[k] * cd_x[k] + g_z_0_x_xyyy[k];

            g_z_0_xx_yyz[k] = -g_z_0_x_yyz[k] * cd_x[k] + g_z_0_x_xyyz[k];

            g_z_0_xx_yzz[k] = -g_z_0_x_yzz[k] * cd_x[k] + g_z_0_x_xyzz[k];

            g_z_0_xx_zzz[k] = -g_z_0_x_zzz[k] * cd_x[k] + g_z_0_x_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_z_0_xy_xxx = cbuffer.data(df_geom_10_off + 120 * acomps  + 10);

        auto g_z_0_xy_xxy = cbuffer.data(df_geom_10_off + 120 * acomps  + 11);

        auto g_z_0_xy_xxz = cbuffer.data(df_geom_10_off + 120 * acomps  + 12);

        auto g_z_0_xy_xyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 13);

        auto g_z_0_xy_xyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 14);

        auto g_z_0_xy_xzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 15);

        auto g_z_0_xy_yyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 16);

        auto g_z_0_xy_yyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 17);

        auto g_z_0_xy_yzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 18);

        auto g_z_0_xy_zzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 19);

        #pragma omp simd aligned(cd_x, g_z_0_xy_xxx, g_z_0_xy_xxy, g_z_0_xy_xxz, g_z_0_xy_xyy, g_z_0_xy_xyz, g_z_0_xy_xzz, g_z_0_xy_yyy, g_z_0_xy_yyz, g_z_0_xy_yzz, g_z_0_xy_zzz, g_z_0_y_xxx, g_z_0_y_xxxx, g_z_0_y_xxxy, g_z_0_y_xxxz, g_z_0_y_xxy, g_z_0_y_xxyy, g_z_0_y_xxyz, g_z_0_y_xxz, g_z_0_y_xxzz, g_z_0_y_xyy, g_z_0_y_xyyy, g_z_0_y_xyyz, g_z_0_y_xyz, g_z_0_y_xyzz, g_z_0_y_xzz, g_z_0_y_xzzz, g_z_0_y_yyy, g_z_0_y_yyz, g_z_0_y_yzz, g_z_0_y_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xy_xxx[k] = -g_z_0_y_xxx[k] * cd_x[k] + g_z_0_y_xxxx[k];

            g_z_0_xy_xxy[k] = -g_z_0_y_xxy[k] * cd_x[k] + g_z_0_y_xxxy[k];

            g_z_0_xy_xxz[k] = -g_z_0_y_xxz[k] * cd_x[k] + g_z_0_y_xxxz[k];

            g_z_0_xy_xyy[k] = -g_z_0_y_xyy[k] * cd_x[k] + g_z_0_y_xxyy[k];

            g_z_0_xy_xyz[k] = -g_z_0_y_xyz[k] * cd_x[k] + g_z_0_y_xxyz[k];

            g_z_0_xy_xzz[k] = -g_z_0_y_xzz[k] * cd_x[k] + g_z_0_y_xxzz[k];

            g_z_0_xy_yyy[k] = -g_z_0_y_yyy[k] * cd_x[k] + g_z_0_y_xyyy[k];

            g_z_0_xy_yyz[k] = -g_z_0_y_yyz[k] * cd_x[k] + g_z_0_y_xyyz[k];

            g_z_0_xy_yzz[k] = -g_z_0_y_yzz[k] * cd_x[k] + g_z_0_y_xyzz[k];

            g_z_0_xy_zzz[k] = -g_z_0_y_zzz[k] * cd_x[k] + g_z_0_y_xzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_z_0_xz_xxx = cbuffer.data(df_geom_10_off + 120 * acomps  + 20);

        auto g_z_0_xz_xxy = cbuffer.data(df_geom_10_off + 120 * acomps  + 21);

        auto g_z_0_xz_xxz = cbuffer.data(df_geom_10_off + 120 * acomps  + 22);

        auto g_z_0_xz_xyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 23);

        auto g_z_0_xz_xyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 24);

        auto g_z_0_xz_xzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 25);

        auto g_z_0_xz_yyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 26);

        auto g_z_0_xz_yyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 27);

        auto g_z_0_xz_yzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 28);

        auto g_z_0_xz_zzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_z_0_xz_xxx, g_z_0_xz_xxy, g_z_0_xz_xxz, g_z_0_xz_xyy, g_z_0_xz_xyz, g_z_0_xz_xzz, g_z_0_xz_yyy, g_z_0_xz_yyz, g_z_0_xz_yzz, g_z_0_xz_zzz, g_z_0_z_xxx, g_z_0_z_xxxx, g_z_0_z_xxxy, g_z_0_z_xxxz, g_z_0_z_xxy, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxz, g_z_0_z_xxzz, g_z_0_z_xyy, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyz, g_z_0_z_xyzz, g_z_0_z_xzz, g_z_0_z_xzzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xz_xxx[k] = -g_z_0_z_xxx[k] * cd_x[k] + g_z_0_z_xxxx[k];

            g_z_0_xz_xxy[k] = -g_z_0_z_xxy[k] * cd_x[k] + g_z_0_z_xxxy[k];

            g_z_0_xz_xxz[k] = -g_z_0_z_xxz[k] * cd_x[k] + g_z_0_z_xxxz[k];

            g_z_0_xz_xyy[k] = -g_z_0_z_xyy[k] * cd_x[k] + g_z_0_z_xxyy[k];

            g_z_0_xz_xyz[k] = -g_z_0_z_xyz[k] * cd_x[k] + g_z_0_z_xxyz[k];

            g_z_0_xz_xzz[k] = -g_z_0_z_xzz[k] * cd_x[k] + g_z_0_z_xxzz[k];

            g_z_0_xz_yyy[k] = -g_z_0_z_yyy[k] * cd_x[k] + g_z_0_z_xyyy[k];

            g_z_0_xz_yyz[k] = -g_z_0_z_yyz[k] * cd_x[k] + g_z_0_z_xyyz[k];

            g_z_0_xz_yzz[k] = -g_z_0_z_yzz[k] * cd_x[k] + g_z_0_z_xyzz[k];

            g_z_0_xz_zzz[k] = -g_z_0_z_zzz[k] * cd_x[k] + g_z_0_z_xzzz[k];
        }

        /// Set up 30-40 components of targeted buffer : cbuffer.data(

        auto g_z_0_yy_xxx = cbuffer.data(df_geom_10_off + 120 * acomps  + 30);

        auto g_z_0_yy_xxy = cbuffer.data(df_geom_10_off + 120 * acomps  + 31);

        auto g_z_0_yy_xxz = cbuffer.data(df_geom_10_off + 120 * acomps  + 32);

        auto g_z_0_yy_xyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 33);

        auto g_z_0_yy_xyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 34);

        auto g_z_0_yy_xzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 35);

        auto g_z_0_yy_yyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 36);

        auto g_z_0_yy_yyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 37);

        auto g_z_0_yy_yzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 38);

        auto g_z_0_yy_zzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 39);

        #pragma omp simd aligned(cd_y, g_z_0_y_xxx, g_z_0_y_xxxy, g_z_0_y_xxy, g_z_0_y_xxyy, g_z_0_y_xxyz, g_z_0_y_xxz, g_z_0_y_xyy, g_z_0_y_xyyy, g_z_0_y_xyyz, g_z_0_y_xyz, g_z_0_y_xyzz, g_z_0_y_xzz, g_z_0_y_yyy, g_z_0_y_yyyy, g_z_0_y_yyyz, g_z_0_y_yyz, g_z_0_y_yyzz, g_z_0_y_yzz, g_z_0_y_yzzz, g_z_0_y_zzz, g_z_0_yy_xxx, g_z_0_yy_xxy, g_z_0_yy_xxz, g_z_0_yy_xyy, g_z_0_yy_xyz, g_z_0_yy_xzz, g_z_0_yy_yyy, g_z_0_yy_yyz, g_z_0_yy_yzz, g_z_0_yy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yy_xxx[k] = -g_z_0_y_xxx[k] * cd_y[k] + g_z_0_y_xxxy[k];

            g_z_0_yy_xxy[k] = -g_z_0_y_xxy[k] * cd_y[k] + g_z_0_y_xxyy[k];

            g_z_0_yy_xxz[k] = -g_z_0_y_xxz[k] * cd_y[k] + g_z_0_y_xxyz[k];

            g_z_0_yy_xyy[k] = -g_z_0_y_xyy[k] * cd_y[k] + g_z_0_y_xyyy[k];

            g_z_0_yy_xyz[k] = -g_z_0_y_xyz[k] * cd_y[k] + g_z_0_y_xyyz[k];

            g_z_0_yy_xzz[k] = -g_z_0_y_xzz[k] * cd_y[k] + g_z_0_y_xyzz[k];

            g_z_0_yy_yyy[k] = -g_z_0_y_yyy[k] * cd_y[k] + g_z_0_y_yyyy[k];

            g_z_0_yy_yyz[k] = -g_z_0_y_yyz[k] * cd_y[k] + g_z_0_y_yyyz[k];

            g_z_0_yy_yzz[k] = -g_z_0_y_yzz[k] * cd_y[k] + g_z_0_y_yyzz[k];

            g_z_0_yy_zzz[k] = -g_z_0_y_zzz[k] * cd_y[k] + g_z_0_y_yzzz[k];
        }

        /// Set up 40-50 components of targeted buffer : cbuffer.data(

        auto g_z_0_yz_xxx = cbuffer.data(df_geom_10_off + 120 * acomps  + 40);

        auto g_z_0_yz_xxy = cbuffer.data(df_geom_10_off + 120 * acomps  + 41);

        auto g_z_0_yz_xxz = cbuffer.data(df_geom_10_off + 120 * acomps  + 42);

        auto g_z_0_yz_xyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 43);

        auto g_z_0_yz_xyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 44);

        auto g_z_0_yz_xzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 45);

        auto g_z_0_yz_yyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 46);

        auto g_z_0_yz_yyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 47);

        auto g_z_0_yz_yzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 48);

        auto g_z_0_yz_zzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 49);

        #pragma omp simd aligned(cd_y, g_z_0_yz_xxx, g_z_0_yz_xxy, g_z_0_yz_xxz, g_z_0_yz_xyy, g_z_0_yz_xyz, g_z_0_yz_xzz, g_z_0_yz_yyy, g_z_0_yz_yyz, g_z_0_yz_yzz, g_z_0_yz_zzz, g_z_0_z_xxx, g_z_0_z_xxxy, g_z_0_z_xxy, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyz, g_z_0_z_xyzz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyz, g_z_0_z_yyzz, g_z_0_z_yzz, g_z_0_z_yzzz, g_z_0_z_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yz_xxx[k] = -g_z_0_z_xxx[k] * cd_y[k] + g_z_0_z_xxxy[k];

            g_z_0_yz_xxy[k] = -g_z_0_z_xxy[k] * cd_y[k] + g_z_0_z_xxyy[k];

            g_z_0_yz_xxz[k] = -g_z_0_z_xxz[k] * cd_y[k] + g_z_0_z_xxyz[k];

            g_z_0_yz_xyy[k] = -g_z_0_z_xyy[k] * cd_y[k] + g_z_0_z_xyyy[k];

            g_z_0_yz_xyz[k] = -g_z_0_z_xyz[k] * cd_y[k] + g_z_0_z_xyyz[k];

            g_z_0_yz_xzz[k] = -g_z_0_z_xzz[k] * cd_y[k] + g_z_0_z_xyzz[k];

            g_z_0_yz_yyy[k] = -g_z_0_z_yyy[k] * cd_y[k] + g_z_0_z_yyyy[k];

            g_z_0_yz_yyz[k] = -g_z_0_z_yyz[k] * cd_y[k] + g_z_0_z_yyyz[k];

            g_z_0_yz_yzz[k] = -g_z_0_z_yzz[k] * cd_y[k] + g_z_0_z_yyzz[k];

            g_z_0_yz_zzz[k] = -g_z_0_z_zzz[k] * cd_y[k] + g_z_0_z_yzzz[k];
        }

        /// Set up 50-60 components of targeted buffer : cbuffer.data(

        auto g_z_0_zz_xxx = cbuffer.data(df_geom_10_off + 120 * acomps  + 50);

        auto g_z_0_zz_xxy = cbuffer.data(df_geom_10_off + 120 * acomps  + 51);

        auto g_z_0_zz_xxz = cbuffer.data(df_geom_10_off + 120 * acomps  + 52);

        auto g_z_0_zz_xyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 53);

        auto g_z_0_zz_xyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 54);

        auto g_z_0_zz_xzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 55);

        auto g_z_0_zz_yyy = cbuffer.data(df_geom_10_off + 120 * acomps  + 56);

        auto g_z_0_zz_yyz = cbuffer.data(df_geom_10_off + 120 * acomps  + 57);

        auto g_z_0_zz_yzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 58);

        auto g_z_0_zz_zzz = cbuffer.data(df_geom_10_off + 120 * acomps  + 59);

        #pragma omp simd aligned(cd_z, g_z_0_z_xxx, g_z_0_z_xxxz, g_z_0_z_xxy, g_z_0_z_xxyz, g_z_0_z_xxz, g_z_0_z_xxzz, g_z_0_z_xyy, g_z_0_z_xyyz, g_z_0_z_xyz, g_z_0_z_xyzz, g_z_0_z_xzz, g_z_0_z_xzzz, g_z_0_z_yyy, g_z_0_z_yyyz, g_z_0_z_yyz, g_z_0_z_yyzz, g_z_0_z_yzz, g_z_0_z_yzzz, g_z_0_z_zzz, g_z_0_z_zzzz, g_z_0_zz_xxx, g_z_0_zz_xxy, g_z_0_zz_xxz, g_z_0_zz_xyy, g_z_0_zz_xyz, g_z_0_zz_xzz, g_z_0_zz_yyy, g_z_0_zz_yyz, g_z_0_zz_yzz, g_z_0_zz_zzz, g_z_xxx, g_z_xxy, g_z_xxz, g_z_xyy, g_z_xyz, g_z_xzz, g_z_yyy, g_z_yyz, g_z_yzz, g_z_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zz_xxx[k] = -g_z_xxx[k] - g_z_0_z_xxx[k] * cd_z[k] + g_z_0_z_xxxz[k];

            g_z_0_zz_xxy[k] = -g_z_xxy[k] - g_z_0_z_xxy[k] * cd_z[k] + g_z_0_z_xxyz[k];

            g_z_0_zz_xxz[k] = -g_z_xxz[k] - g_z_0_z_xxz[k] * cd_z[k] + g_z_0_z_xxzz[k];

            g_z_0_zz_xyy[k] = -g_z_xyy[k] - g_z_0_z_xyy[k] * cd_z[k] + g_z_0_z_xyyz[k];

            g_z_0_zz_xyz[k] = -g_z_xyz[k] - g_z_0_z_xyz[k] * cd_z[k] + g_z_0_z_xyzz[k];

            g_z_0_zz_xzz[k] = -g_z_xzz[k] - g_z_0_z_xzz[k] * cd_z[k] + g_z_0_z_xzzz[k];

            g_z_0_zz_yyy[k] = -g_z_yyy[k] - g_z_0_z_yyy[k] * cd_z[k] + g_z_0_z_yyyz[k];

            g_z_0_zz_yyz[k] = -g_z_yyz[k] - g_z_0_z_yyz[k] * cd_z[k] + g_z_0_z_yyzz[k];

            g_z_0_zz_yzz[k] = -g_z_yzz[k] - g_z_0_z_yzz[k] * cd_z[k] + g_z_0_z_yzzz[k];

            g_z_0_zz_zzz[k] = -g_z_zzz[k] - g_z_0_z_zzz[k] * cd_z[k] + g_z_0_z_zzzz[k];
        }
    }
}

} // t3ceri namespace

