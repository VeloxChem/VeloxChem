#include "ElectronRepulsionGeom2000ContrRecPFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_pfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_pfxx,
                                            const size_t idx_geom_10_sfxx,
                                            const size_t idx_geom_20_sfxx,
                                            const size_t idx_geom_20_sgxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SFSS

            const auto sf_geom_10_off = idx_geom_10_sfxx + i * dcomps + j;

            auto g_x_0_0_xxx = cbuffer.data(sf_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxy = cbuffer.data(sf_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxz = cbuffer.data(sf_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xyy = cbuffer.data(sf_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xyz = cbuffer.data(sf_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xzz = cbuffer.data(sf_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_0_yyy = cbuffer.data(sf_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_0_yyz = cbuffer.data(sf_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_0_yzz = cbuffer.data(sf_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_0_zzz = cbuffer.data(sf_geom_10_off + 9 * ccomps * dcomps);

            auto g_y_0_0_xxx = cbuffer.data(sf_geom_10_off + 10 * ccomps * dcomps);

            auto g_y_0_0_xxy = cbuffer.data(sf_geom_10_off + 11 * ccomps * dcomps);

            auto g_y_0_0_xxz = cbuffer.data(sf_geom_10_off + 12 * ccomps * dcomps);

            auto g_y_0_0_xyy = cbuffer.data(sf_geom_10_off + 13 * ccomps * dcomps);

            auto g_y_0_0_xyz = cbuffer.data(sf_geom_10_off + 14 * ccomps * dcomps);

            auto g_y_0_0_xzz = cbuffer.data(sf_geom_10_off + 15 * ccomps * dcomps);

            auto g_y_0_0_yyy = cbuffer.data(sf_geom_10_off + 16 * ccomps * dcomps);

            auto g_y_0_0_yyz = cbuffer.data(sf_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_0_yzz = cbuffer.data(sf_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_0_zzz = cbuffer.data(sf_geom_10_off + 19 * ccomps * dcomps);

            auto g_z_0_0_xxx = cbuffer.data(sf_geom_10_off + 20 * ccomps * dcomps);

            auto g_z_0_0_xxy = cbuffer.data(sf_geom_10_off + 21 * ccomps * dcomps);

            auto g_z_0_0_xxz = cbuffer.data(sf_geom_10_off + 22 * ccomps * dcomps);

            auto g_z_0_0_xyy = cbuffer.data(sf_geom_10_off + 23 * ccomps * dcomps);

            auto g_z_0_0_xyz = cbuffer.data(sf_geom_10_off + 24 * ccomps * dcomps);

            auto g_z_0_0_xzz = cbuffer.data(sf_geom_10_off + 25 * ccomps * dcomps);

            auto g_z_0_0_yyy = cbuffer.data(sf_geom_10_off + 26 * ccomps * dcomps);

            auto g_z_0_0_yyz = cbuffer.data(sf_geom_10_off + 27 * ccomps * dcomps);

            auto g_z_0_0_yzz = cbuffer.data(sf_geom_10_off + 28 * ccomps * dcomps);

            auto g_z_0_0_zzz = cbuffer.data(sf_geom_10_off + 29 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SFSS

            const auto sf_geom_20_off = idx_geom_20_sfxx + i * dcomps + j;

            auto g_xx_0_0_xxx = cbuffer.data(sf_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_0_xxy = cbuffer.data(sf_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_0_xxz = cbuffer.data(sf_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_0_xyy = cbuffer.data(sf_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_0_xyz = cbuffer.data(sf_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_0_xzz = cbuffer.data(sf_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_0_yyy = cbuffer.data(sf_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_0_yyz = cbuffer.data(sf_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_0_yzz = cbuffer.data(sf_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_0_zzz = cbuffer.data(sf_geom_20_off + 9 * ccomps * dcomps);

            auto g_xy_0_0_xxx = cbuffer.data(sf_geom_20_off + 10 * ccomps * dcomps);

            auto g_xy_0_0_xxy = cbuffer.data(sf_geom_20_off + 11 * ccomps * dcomps);

            auto g_xy_0_0_xxz = cbuffer.data(sf_geom_20_off + 12 * ccomps * dcomps);

            auto g_xy_0_0_xyy = cbuffer.data(sf_geom_20_off + 13 * ccomps * dcomps);

            auto g_xy_0_0_xyz = cbuffer.data(sf_geom_20_off + 14 * ccomps * dcomps);

            auto g_xy_0_0_xzz = cbuffer.data(sf_geom_20_off + 15 * ccomps * dcomps);

            auto g_xy_0_0_yyy = cbuffer.data(sf_geom_20_off + 16 * ccomps * dcomps);

            auto g_xy_0_0_yyz = cbuffer.data(sf_geom_20_off + 17 * ccomps * dcomps);

            auto g_xy_0_0_yzz = cbuffer.data(sf_geom_20_off + 18 * ccomps * dcomps);

            auto g_xy_0_0_zzz = cbuffer.data(sf_geom_20_off + 19 * ccomps * dcomps);

            auto g_xz_0_0_xxx = cbuffer.data(sf_geom_20_off + 20 * ccomps * dcomps);

            auto g_xz_0_0_xxy = cbuffer.data(sf_geom_20_off + 21 * ccomps * dcomps);

            auto g_xz_0_0_xxz = cbuffer.data(sf_geom_20_off + 22 * ccomps * dcomps);

            auto g_xz_0_0_xyy = cbuffer.data(sf_geom_20_off + 23 * ccomps * dcomps);

            auto g_xz_0_0_xyz = cbuffer.data(sf_geom_20_off + 24 * ccomps * dcomps);

            auto g_xz_0_0_xzz = cbuffer.data(sf_geom_20_off + 25 * ccomps * dcomps);

            auto g_xz_0_0_yyy = cbuffer.data(sf_geom_20_off + 26 * ccomps * dcomps);

            auto g_xz_0_0_yyz = cbuffer.data(sf_geom_20_off + 27 * ccomps * dcomps);

            auto g_xz_0_0_yzz = cbuffer.data(sf_geom_20_off + 28 * ccomps * dcomps);

            auto g_xz_0_0_zzz = cbuffer.data(sf_geom_20_off + 29 * ccomps * dcomps);

            auto g_yy_0_0_xxx = cbuffer.data(sf_geom_20_off + 30 * ccomps * dcomps);

            auto g_yy_0_0_xxy = cbuffer.data(sf_geom_20_off + 31 * ccomps * dcomps);

            auto g_yy_0_0_xxz = cbuffer.data(sf_geom_20_off + 32 * ccomps * dcomps);

            auto g_yy_0_0_xyy = cbuffer.data(sf_geom_20_off + 33 * ccomps * dcomps);

            auto g_yy_0_0_xyz = cbuffer.data(sf_geom_20_off + 34 * ccomps * dcomps);

            auto g_yy_0_0_xzz = cbuffer.data(sf_geom_20_off + 35 * ccomps * dcomps);

            auto g_yy_0_0_yyy = cbuffer.data(sf_geom_20_off + 36 * ccomps * dcomps);

            auto g_yy_0_0_yyz = cbuffer.data(sf_geom_20_off + 37 * ccomps * dcomps);

            auto g_yy_0_0_yzz = cbuffer.data(sf_geom_20_off + 38 * ccomps * dcomps);

            auto g_yy_0_0_zzz = cbuffer.data(sf_geom_20_off + 39 * ccomps * dcomps);

            auto g_yz_0_0_xxx = cbuffer.data(sf_geom_20_off + 40 * ccomps * dcomps);

            auto g_yz_0_0_xxy = cbuffer.data(sf_geom_20_off + 41 * ccomps * dcomps);

            auto g_yz_0_0_xxz = cbuffer.data(sf_geom_20_off + 42 * ccomps * dcomps);

            auto g_yz_0_0_xyy = cbuffer.data(sf_geom_20_off + 43 * ccomps * dcomps);

            auto g_yz_0_0_xyz = cbuffer.data(sf_geom_20_off + 44 * ccomps * dcomps);

            auto g_yz_0_0_xzz = cbuffer.data(sf_geom_20_off + 45 * ccomps * dcomps);

            auto g_yz_0_0_yyy = cbuffer.data(sf_geom_20_off + 46 * ccomps * dcomps);

            auto g_yz_0_0_yyz = cbuffer.data(sf_geom_20_off + 47 * ccomps * dcomps);

            auto g_yz_0_0_yzz = cbuffer.data(sf_geom_20_off + 48 * ccomps * dcomps);

            auto g_yz_0_0_zzz = cbuffer.data(sf_geom_20_off + 49 * ccomps * dcomps);

            auto g_zz_0_0_xxx = cbuffer.data(sf_geom_20_off + 50 * ccomps * dcomps);

            auto g_zz_0_0_xxy = cbuffer.data(sf_geom_20_off + 51 * ccomps * dcomps);

            auto g_zz_0_0_xxz = cbuffer.data(sf_geom_20_off + 52 * ccomps * dcomps);

            auto g_zz_0_0_xyy = cbuffer.data(sf_geom_20_off + 53 * ccomps * dcomps);

            auto g_zz_0_0_xyz = cbuffer.data(sf_geom_20_off + 54 * ccomps * dcomps);

            auto g_zz_0_0_xzz = cbuffer.data(sf_geom_20_off + 55 * ccomps * dcomps);

            auto g_zz_0_0_yyy = cbuffer.data(sf_geom_20_off + 56 * ccomps * dcomps);

            auto g_zz_0_0_yyz = cbuffer.data(sf_geom_20_off + 57 * ccomps * dcomps);

            auto g_zz_0_0_yzz = cbuffer.data(sf_geom_20_off + 58 * ccomps * dcomps);

            auto g_zz_0_0_zzz = cbuffer.data(sf_geom_20_off + 59 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SGSS

            const auto sg_geom_20_off = idx_geom_20_sgxx + i * dcomps + j;

            auto g_xx_0_0_xxxx = cbuffer.data(sg_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_0_xxxy = cbuffer.data(sg_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_0_xxxz = cbuffer.data(sg_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_0_xxyy = cbuffer.data(sg_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_0_xxyz = cbuffer.data(sg_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_0_xxzz = cbuffer.data(sg_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_0_xyyy = cbuffer.data(sg_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_0_xyyz = cbuffer.data(sg_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_0_xyzz = cbuffer.data(sg_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_0_xzzz = cbuffer.data(sg_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_0_yyyy = cbuffer.data(sg_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_0_yyyz = cbuffer.data(sg_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_0_yyzz = cbuffer.data(sg_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_0_yzzz = cbuffer.data(sg_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_0_zzzz = cbuffer.data(sg_geom_20_off + 14 * ccomps * dcomps);

            auto g_xy_0_0_xxxx = cbuffer.data(sg_geom_20_off + 15 * ccomps * dcomps);

            auto g_xy_0_0_xxxy = cbuffer.data(sg_geom_20_off + 16 * ccomps * dcomps);

            auto g_xy_0_0_xxxz = cbuffer.data(sg_geom_20_off + 17 * ccomps * dcomps);

            auto g_xy_0_0_xxyy = cbuffer.data(sg_geom_20_off + 18 * ccomps * dcomps);

            auto g_xy_0_0_xxyz = cbuffer.data(sg_geom_20_off + 19 * ccomps * dcomps);

            auto g_xy_0_0_xxzz = cbuffer.data(sg_geom_20_off + 20 * ccomps * dcomps);

            auto g_xy_0_0_xyyy = cbuffer.data(sg_geom_20_off + 21 * ccomps * dcomps);

            auto g_xy_0_0_xyyz = cbuffer.data(sg_geom_20_off + 22 * ccomps * dcomps);

            auto g_xy_0_0_xyzz = cbuffer.data(sg_geom_20_off + 23 * ccomps * dcomps);

            auto g_xy_0_0_xzzz = cbuffer.data(sg_geom_20_off + 24 * ccomps * dcomps);

            auto g_xy_0_0_yyyy = cbuffer.data(sg_geom_20_off + 25 * ccomps * dcomps);

            auto g_xy_0_0_yyyz = cbuffer.data(sg_geom_20_off + 26 * ccomps * dcomps);

            auto g_xy_0_0_yyzz = cbuffer.data(sg_geom_20_off + 27 * ccomps * dcomps);

            auto g_xy_0_0_yzzz = cbuffer.data(sg_geom_20_off + 28 * ccomps * dcomps);

            auto g_xy_0_0_zzzz = cbuffer.data(sg_geom_20_off + 29 * ccomps * dcomps);

            auto g_xz_0_0_xxxx = cbuffer.data(sg_geom_20_off + 30 * ccomps * dcomps);

            auto g_xz_0_0_xxxy = cbuffer.data(sg_geom_20_off + 31 * ccomps * dcomps);

            auto g_xz_0_0_xxxz = cbuffer.data(sg_geom_20_off + 32 * ccomps * dcomps);

            auto g_xz_0_0_xxyy = cbuffer.data(sg_geom_20_off + 33 * ccomps * dcomps);

            auto g_xz_0_0_xxyz = cbuffer.data(sg_geom_20_off + 34 * ccomps * dcomps);

            auto g_xz_0_0_xxzz = cbuffer.data(sg_geom_20_off + 35 * ccomps * dcomps);

            auto g_xz_0_0_xyyy = cbuffer.data(sg_geom_20_off + 36 * ccomps * dcomps);

            auto g_xz_0_0_xyyz = cbuffer.data(sg_geom_20_off + 37 * ccomps * dcomps);

            auto g_xz_0_0_xyzz = cbuffer.data(sg_geom_20_off + 38 * ccomps * dcomps);

            auto g_xz_0_0_xzzz = cbuffer.data(sg_geom_20_off + 39 * ccomps * dcomps);

            auto g_xz_0_0_yyyy = cbuffer.data(sg_geom_20_off + 40 * ccomps * dcomps);

            auto g_xz_0_0_yyyz = cbuffer.data(sg_geom_20_off + 41 * ccomps * dcomps);

            auto g_xz_0_0_yyzz = cbuffer.data(sg_geom_20_off + 42 * ccomps * dcomps);

            auto g_xz_0_0_yzzz = cbuffer.data(sg_geom_20_off + 43 * ccomps * dcomps);

            auto g_xz_0_0_zzzz = cbuffer.data(sg_geom_20_off + 44 * ccomps * dcomps);

            auto g_yy_0_0_xxxx = cbuffer.data(sg_geom_20_off + 45 * ccomps * dcomps);

            auto g_yy_0_0_xxxy = cbuffer.data(sg_geom_20_off + 46 * ccomps * dcomps);

            auto g_yy_0_0_xxxz = cbuffer.data(sg_geom_20_off + 47 * ccomps * dcomps);

            auto g_yy_0_0_xxyy = cbuffer.data(sg_geom_20_off + 48 * ccomps * dcomps);

            auto g_yy_0_0_xxyz = cbuffer.data(sg_geom_20_off + 49 * ccomps * dcomps);

            auto g_yy_0_0_xxzz = cbuffer.data(sg_geom_20_off + 50 * ccomps * dcomps);

            auto g_yy_0_0_xyyy = cbuffer.data(sg_geom_20_off + 51 * ccomps * dcomps);

            auto g_yy_0_0_xyyz = cbuffer.data(sg_geom_20_off + 52 * ccomps * dcomps);

            auto g_yy_0_0_xyzz = cbuffer.data(sg_geom_20_off + 53 * ccomps * dcomps);

            auto g_yy_0_0_xzzz = cbuffer.data(sg_geom_20_off + 54 * ccomps * dcomps);

            auto g_yy_0_0_yyyy = cbuffer.data(sg_geom_20_off + 55 * ccomps * dcomps);

            auto g_yy_0_0_yyyz = cbuffer.data(sg_geom_20_off + 56 * ccomps * dcomps);

            auto g_yy_0_0_yyzz = cbuffer.data(sg_geom_20_off + 57 * ccomps * dcomps);

            auto g_yy_0_0_yzzz = cbuffer.data(sg_geom_20_off + 58 * ccomps * dcomps);

            auto g_yy_0_0_zzzz = cbuffer.data(sg_geom_20_off + 59 * ccomps * dcomps);

            auto g_yz_0_0_xxxx = cbuffer.data(sg_geom_20_off + 60 * ccomps * dcomps);

            auto g_yz_0_0_xxxy = cbuffer.data(sg_geom_20_off + 61 * ccomps * dcomps);

            auto g_yz_0_0_xxxz = cbuffer.data(sg_geom_20_off + 62 * ccomps * dcomps);

            auto g_yz_0_0_xxyy = cbuffer.data(sg_geom_20_off + 63 * ccomps * dcomps);

            auto g_yz_0_0_xxyz = cbuffer.data(sg_geom_20_off + 64 * ccomps * dcomps);

            auto g_yz_0_0_xxzz = cbuffer.data(sg_geom_20_off + 65 * ccomps * dcomps);

            auto g_yz_0_0_xyyy = cbuffer.data(sg_geom_20_off + 66 * ccomps * dcomps);

            auto g_yz_0_0_xyyz = cbuffer.data(sg_geom_20_off + 67 * ccomps * dcomps);

            auto g_yz_0_0_xyzz = cbuffer.data(sg_geom_20_off + 68 * ccomps * dcomps);

            auto g_yz_0_0_xzzz = cbuffer.data(sg_geom_20_off + 69 * ccomps * dcomps);

            auto g_yz_0_0_yyyy = cbuffer.data(sg_geom_20_off + 70 * ccomps * dcomps);

            auto g_yz_0_0_yyyz = cbuffer.data(sg_geom_20_off + 71 * ccomps * dcomps);

            auto g_yz_0_0_yyzz = cbuffer.data(sg_geom_20_off + 72 * ccomps * dcomps);

            auto g_yz_0_0_yzzz = cbuffer.data(sg_geom_20_off + 73 * ccomps * dcomps);

            auto g_yz_0_0_zzzz = cbuffer.data(sg_geom_20_off + 74 * ccomps * dcomps);

            auto g_zz_0_0_xxxx = cbuffer.data(sg_geom_20_off + 75 * ccomps * dcomps);

            auto g_zz_0_0_xxxy = cbuffer.data(sg_geom_20_off + 76 * ccomps * dcomps);

            auto g_zz_0_0_xxxz = cbuffer.data(sg_geom_20_off + 77 * ccomps * dcomps);

            auto g_zz_0_0_xxyy = cbuffer.data(sg_geom_20_off + 78 * ccomps * dcomps);

            auto g_zz_0_0_xxyz = cbuffer.data(sg_geom_20_off + 79 * ccomps * dcomps);

            auto g_zz_0_0_xxzz = cbuffer.data(sg_geom_20_off + 80 * ccomps * dcomps);

            auto g_zz_0_0_xyyy = cbuffer.data(sg_geom_20_off + 81 * ccomps * dcomps);

            auto g_zz_0_0_xyyz = cbuffer.data(sg_geom_20_off + 82 * ccomps * dcomps);

            auto g_zz_0_0_xyzz = cbuffer.data(sg_geom_20_off + 83 * ccomps * dcomps);

            auto g_zz_0_0_xzzz = cbuffer.data(sg_geom_20_off + 84 * ccomps * dcomps);

            auto g_zz_0_0_yyyy = cbuffer.data(sg_geom_20_off + 85 * ccomps * dcomps);

            auto g_zz_0_0_yyyz = cbuffer.data(sg_geom_20_off + 86 * ccomps * dcomps);

            auto g_zz_0_0_yyzz = cbuffer.data(sg_geom_20_off + 87 * ccomps * dcomps);

            auto g_zz_0_0_yzzz = cbuffer.data(sg_geom_20_off + 88 * ccomps * dcomps);

            auto g_zz_0_0_zzzz = cbuffer.data(sg_geom_20_off + 89 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pfxx

            const auto pf_geom_20_off = idx_geom_20_pfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_xx_0_x_xxx = cbuffer.data(pf_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_x_xxy = cbuffer.data(pf_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_x_xxz = cbuffer.data(pf_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_x_xyy = cbuffer.data(pf_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_x_xyz = cbuffer.data(pf_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_x_xzz = cbuffer.data(pf_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_x_yyy = cbuffer.data(pf_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_x_yyz = cbuffer.data(pf_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_x_yzz = cbuffer.data(pf_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_x_zzz = cbuffer.data(pf_geom_20_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxx, g_x_0_0_xxy, g_x_0_0_xxz, g_x_0_0_xyy, g_x_0_0_xyz, g_x_0_0_xzz, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yzz, g_x_0_0_zzz, g_xx_0_0_xxx, g_xx_0_0_xxxx, g_xx_0_0_xxxy, g_xx_0_0_xxxz, g_xx_0_0_xxy, g_xx_0_0_xxyy, g_xx_0_0_xxyz, g_xx_0_0_xxz, g_xx_0_0_xxzz, g_xx_0_0_xyy, g_xx_0_0_xyyy, g_xx_0_0_xyyz, g_xx_0_0_xyz, g_xx_0_0_xyzz, g_xx_0_0_xzz, g_xx_0_0_xzzz, g_xx_0_0_yyy, g_xx_0_0_yyz, g_xx_0_0_yzz, g_xx_0_0_zzz, g_xx_0_x_xxx, g_xx_0_x_xxy, g_xx_0_x_xxz, g_xx_0_x_xyy, g_xx_0_x_xyz, g_xx_0_x_xzz, g_xx_0_x_yyy, g_xx_0_x_yyz, g_xx_0_x_yzz, g_xx_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_x_xxx[k] = -2.0 * g_x_0_0_xxx[k] - g_xx_0_0_xxx[k] * ab_x + g_xx_0_0_xxxx[k];

                g_xx_0_x_xxy[k] = -2.0 * g_x_0_0_xxy[k] - g_xx_0_0_xxy[k] * ab_x + g_xx_0_0_xxxy[k];

                g_xx_0_x_xxz[k] = -2.0 * g_x_0_0_xxz[k] - g_xx_0_0_xxz[k] * ab_x + g_xx_0_0_xxxz[k];

                g_xx_0_x_xyy[k] = -2.0 * g_x_0_0_xyy[k] - g_xx_0_0_xyy[k] * ab_x + g_xx_0_0_xxyy[k];

                g_xx_0_x_xyz[k] = -2.0 * g_x_0_0_xyz[k] - g_xx_0_0_xyz[k] * ab_x + g_xx_0_0_xxyz[k];

                g_xx_0_x_xzz[k] = -2.0 * g_x_0_0_xzz[k] - g_xx_0_0_xzz[k] * ab_x + g_xx_0_0_xxzz[k];

                g_xx_0_x_yyy[k] = -2.0 * g_x_0_0_yyy[k] - g_xx_0_0_yyy[k] * ab_x + g_xx_0_0_xyyy[k];

                g_xx_0_x_yyz[k] = -2.0 * g_x_0_0_yyz[k] - g_xx_0_0_yyz[k] * ab_x + g_xx_0_0_xyyz[k];

                g_xx_0_x_yzz[k] = -2.0 * g_x_0_0_yzz[k] - g_xx_0_0_yzz[k] * ab_x + g_xx_0_0_xyzz[k];

                g_xx_0_x_zzz[k] = -2.0 * g_x_0_0_zzz[k] - g_xx_0_0_zzz[k] * ab_x + g_xx_0_0_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_xx_0_y_xxx = cbuffer.data(pf_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_y_xxy = cbuffer.data(pf_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_y_xxz = cbuffer.data(pf_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_y_xyy = cbuffer.data(pf_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_y_xyz = cbuffer.data(pf_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_y_xzz = cbuffer.data(pf_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_y_yyy = cbuffer.data(pf_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_y_yyz = cbuffer.data(pf_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_y_yzz = cbuffer.data(pf_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_y_zzz = cbuffer.data(pf_geom_20_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_0_xxx, g_xx_0_0_xxxy, g_xx_0_0_xxy, g_xx_0_0_xxyy, g_xx_0_0_xxyz, g_xx_0_0_xxz, g_xx_0_0_xyy, g_xx_0_0_xyyy, g_xx_0_0_xyyz, g_xx_0_0_xyz, g_xx_0_0_xyzz, g_xx_0_0_xzz, g_xx_0_0_yyy, g_xx_0_0_yyyy, g_xx_0_0_yyyz, g_xx_0_0_yyz, g_xx_0_0_yyzz, g_xx_0_0_yzz, g_xx_0_0_yzzz, g_xx_0_0_zzz, g_xx_0_y_xxx, g_xx_0_y_xxy, g_xx_0_y_xxz, g_xx_0_y_xyy, g_xx_0_y_xyz, g_xx_0_y_xzz, g_xx_0_y_yyy, g_xx_0_y_yyz, g_xx_0_y_yzz, g_xx_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_y_xxx[k] = -g_xx_0_0_xxx[k] * ab_y + g_xx_0_0_xxxy[k];

                g_xx_0_y_xxy[k] = -g_xx_0_0_xxy[k] * ab_y + g_xx_0_0_xxyy[k];

                g_xx_0_y_xxz[k] = -g_xx_0_0_xxz[k] * ab_y + g_xx_0_0_xxyz[k];

                g_xx_0_y_xyy[k] = -g_xx_0_0_xyy[k] * ab_y + g_xx_0_0_xyyy[k];

                g_xx_0_y_xyz[k] = -g_xx_0_0_xyz[k] * ab_y + g_xx_0_0_xyyz[k];

                g_xx_0_y_xzz[k] = -g_xx_0_0_xzz[k] * ab_y + g_xx_0_0_xyzz[k];

                g_xx_0_y_yyy[k] = -g_xx_0_0_yyy[k] * ab_y + g_xx_0_0_yyyy[k];

                g_xx_0_y_yyz[k] = -g_xx_0_0_yyz[k] * ab_y + g_xx_0_0_yyyz[k];

                g_xx_0_y_yzz[k] = -g_xx_0_0_yzz[k] * ab_y + g_xx_0_0_yyzz[k];

                g_xx_0_y_zzz[k] = -g_xx_0_0_zzz[k] * ab_y + g_xx_0_0_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_z_xxx = cbuffer.data(pf_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_z_xxy = cbuffer.data(pf_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_z_xxz = cbuffer.data(pf_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_z_xyy = cbuffer.data(pf_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_z_xyz = cbuffer.data(pf_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_z_xzz = cbuffer.data(pf_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_z_yyy = cbuffer.data(pf_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_z_yyz = cbuffer.data(pf_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_z_yzz = cbuffer.data(pf_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_z_zzz = cbuffer.data(pf_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_0_xxx, g_xx_0_0_xxxz, g_xx_0_0_xxy, g_xx_0_0_xxyz, g_xx_0_0_xxz, g_xx_0_0_xxzz, g_xx_0_0_xyy, g_xx_0_0_xyyz, g_xx_0_0_xyz, g_xx_0_0_xyzz, g_xx_0_0_xzz, g_xx_0_0_xzzz, g_xx_0_0_yyy, g_xx_0_0_yyyz, g_xx_0_0_yyz, g_xx_0_0_yyzz, g_xx_0_0_yzz, g_xx_0_0_yzzz, g_xx_0_0_zzz, g_xx_0_0_zzzz, g_xx_0_z_xxx, g_xx_0_z_xxy, g_xx_0_z_xxz, g_xx_0_z_xyy, g_xx_0_z_xyz, g_xx_0_z_xzz, g_xx_0_z_yyy, g_xx_0_z_yyz, g_xx_0_z_yzz, g_xx_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_z_xxx[k] = -g_xx_0_0_xxx[k] * ab_z + g_xx_0_0_xxxz[k];

                g_xx_0_z_xxy[k] = -g_xx_0_0_xxy[k] * ab_z + g_xx_0_0_xxyz[k];

                g_xx_0_z_xxz[k] = -g_xx_0_0_xxz[k] * ab_z + g_xx_0_0_xxzz[k];

                g_xx_0_z_xyy[k] = -g_xx_0_0_xyy[k] * ab_z + g_xx_0_0_xyyz[k];

                g_xx_0_z_xyz[k] = -g_xx_0_0_xyz[k] * ab_z + g_xx_0_0_xyzz[k];

                g_xx_0_z_xzz[k] = -g_xx_0_0_xzz[k] * ab_z + g_xx_0_0_xzzz[k];

                g_xx_0_z_yyy[k] = -g_xx_0_0_yyy[k] * ab_z + g_xx_0_0_yyyz[k];

                g_xx_0_z_yyz[k] = -g_xx_0_0_yyz[k] * ab_z + g_xx_0_0_yyzz[k];

                g_xx_0_z_yzz[k] = -g_xx_0_0_yzz[k] * ab_z + g_xx_0_0_yzzz[k];

                g_xx_0_z_zzz[k] = -g_xx_0_0_zzz[k] * ab_z + g_xx_0_0_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_xy_0_x_xxx = cbuffer.data(pf_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_x_xxy = cbuffer.data(pf_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_x_xxz = cbuffer.data(pf_geom_20_off + 32 * ccomps * dcomps);

            auto g_xy_0_x_xyy = cbuffer.data(pf_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_x_xyz = cbuffer.data(pf_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_x_xzz = cbuffer.data(pf_geom_20_off + 35 * ccomps * dcomps);

            auto g_xy_0_x_yyy = cbuffer.data(pf_geom_20_off + 36 * ccomps * dcomps);

            auto g_xy_0_x_yyz = cbuffer.data(pf_geom_20_off + 37 * ccomps * dcomps);

            auto g_xy_0_x_yzz = cbuffer.data(pf_geom_20_off + 38 * ccomps * dcomps);

            auto g_xy_0_x_zzz = cbuffer.data(pf_geom_20_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_0_xxx, g_xy_0_0_xxxx, g_xy_0_0_xxxy, g_xy_0_0_xxxz, g_xy_0_0_xxy, g_xy_0_0_xxyy, g_xy_0_0_xxyz, g_xy_0_0_xxz, g_xy_0_0_xxzz, g_xy_0_0_xyy, g_xy_0_0_xyyy, g_xy_0_0_xyyz, g_xy_0_0_xyz, g_xy_0_0_xyzz, g_xy_0_0_xzz, g_xy_0_0_xzzz, g_xy_0_0_yyy, g_xy_0_0_yyz, g_xy_0_0_yzz, g_xy_0_0_zzz, g_xy_0_x_xxx, g_xy_0_x_xxy, g_xy_0_x_xxz, g_xy_0_x_xyy, g_xy_0_x_xyz, g_xy_0_x_xzz, g_xy_0_x_yyy, g_xy_0_x_yyz, g_xy_0_x_yzz, g_xy_0_x_zzz, g_y_0_0_xxx, g_y_0_0_xxy, g_y_0_0_xxz, g_y_0_0_xyy, g_y_0_0_xyz, g_y_0_0_xzz, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yzz, g_y_0_0_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_x_xxx[k] = -g_y_0_0_xxx[k] - g_xy_0_0_xxx[k] * ab_x + g_xy_0_0_xxxx[k];

                g_xy_0_x_xxy[k] = -g_y_0_0_xxy[k] - g_xy_0_0_xxy[k] * ab_x + g_xy_0_0_xxxy[k];

                g_xy_0_x_xxz[k] = -g_y_0_0_xxz[k] - g_xy_0_0_xxz[k] * ab_x + g_xy_0_0_xxxz[k];

                g_xy_0_x_xyy[k] = -g_y_0_0_xyy[k] - g_xy_0_0_xyy[k] * ab_x + g_xy_0_0_xxyy[k];

                g_xy_0_x_xyz[k] = -g_y_0_0_xyz[k] - g_xy_0_0_xyz[k] * ab_x + g_xy_0_0_xxyz[k];

                g_xy_0_x_xzz[k] = -g_y_0_0_xzz[k] - g_xy_0_0_xzz[k] * ab_x + g_xy_0_0_xxzz[k];

                g_xy_0_x_yyy[k] = -g_y_0_0_yyy[k] - g_xy_0_0_yyy[k] * ab_x + g_xy_0_0_xyyy[k];

                g_xy_0_x_yyz[k] = -g_y_0_0_yyz[k] - g_xy_0_0_yyz[k] * ab_x + g_xy_0_0_xyyz[k];

                g_xy_0_x_yzz[k] = -g_y_0_0_yzz[k] - g_xy_0_0_yzz[k] * ab_x + g_xy_0_0_xyzz[k];

                g_xy_0_x_zzz[k] = -g_y_0_0_zzz[k] - g_xy_0_0_zzz[k] * ab_x + g_xy_0_0_xzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_xy_0_y_xxx = cbuffer.data(pf_geom_20_off + 40 * ccomps * dcomps);

            auto g_xy_0_y_xxy = cbuffer.data(pf_geom_20_off + 41 * ccomps * dcomps);

            auto g_xy_0_y_xxz = cbuffer.data(pf_geom_20_off + 42 * ccomps * dcomps);

            auto g_xy_0_y_xyy = cbuffer.data(pf_geom_20_off + 43 * ccomps * dcomps);

            auto g_xy_0_y_xyz = cbuffer.data(pf_geom_20_off + 44 * ccomps * dcomps);

            auto g_xy_0_y_xzz = cbuffer.data(pf_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_y_yyy = cbuffer.data(pf_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_y_yyz = cbuffer.data(pf_geom_20_off + 47 * ccomps * dcomps);

            auto g_xy_0_y_yzz = cbuffer.data(pf_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_y_zzz = cbuffer.data(pf_geom_20_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxx, g_x_0_0_xxy, g_x_0_0_xxz, g_x_0_0_xyy, g_x_0_0_xyz, g_x_0_0_xzz, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yzz, g_x_0_0_zzz, g_xy_0_0_xxx, g_xy_0_0_xxxy, g_xy_0_0_xxy, g_xy_0_0_xxyy, g_xy_0_0_xxyz, g_xy_0_0_xxz, g_xy_0_0_xyy, g_xy_0_0_xyyy, g_xy_0_0_xyyz, g_xy_0_0_xyz, g_xy_0_0_xyzz, g_xy_0_0_xzz, g_xy_0_0_yyy, g_xy_0_0_yyyy, g_xy_0_0_yyyz, g_xy_0_0_yyz, g_xy_0_0_yyzz, g_xy_0_0_yzz, g_xy_0_0_yzzz, g_xy_0_0_zzz, g_xy_0_y_xxx, g_xy_0_y_xxy, g_xy_0_y_xxz, g_xy_0_y_xyy, g_xy_0_y_xyz, g_xy_0_y_xzz, g_xy_0_y_yyy, g_xy_0_y_yyz, g_xy_0_y_yzz, g_xy_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_y_xxx[k] = -g_x_0_0_xxx[k] - g_xy_0_0_xxx[k] * ab_y + g_xy_0_0_xxxy[k];

                g_xy_0_y_xxy[k] = -g_x_0_0_xxy[k] - g_xy_0_0_xxy[k] * ab_y + g_xy_0_0_xxyy[k];

                g_xy_0_y_xxz[k] = -g_x_0_0_xxz[k] - g_xy_0_0_xxz[k] * ab_y + g_xy_0_0_xxyz[k];

                g_xy_0_y_xyy[k] = -g_x_0_0_xyy[k] - g_xy_0_0_xyy[k] * ab_y + g_xy_0_0_xyyy[k];

                g_xy_0_y_xyz[k] = -g_x_0_0_xyz[k] - g_xy_0_0_xyz[k] * ab_y + g_xy_0_0_xyyz[k];

                g_xy_0_y_xzz[k] = -g_x_0_0_xzz[k] - g_xy_0_0_xzz[k] * ab_y + g_xy_0_0_xyzz[k];

                g_xy_0_y_yyy[k] = -g_x_0_0_yyy[k] - g_xy_0_0_yyy[k] * ab_y + g_xy_0_0_yyyy[k];

                g_xy_0_y_yyz[k] = -g_x_0_0_yyz[k] - g_xy_0_0_yyz[k] * ab_y + g_xy_0_0_yyyz[k];

                g_xy_0_y_yzz[k] = -g_x_0_0_yzz[k] - g_xy_0_0_yzz[k] * ab_y + g_xy_0_0_yyzz[k];

                g_xy_0_y_zzz[k] = -g_x_0_0_zzz[k] - g_xy_0_0_zzz[k] * ab_y + g_xy_0_0_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_xy_0_z_xxx = cbuffer.data(pf_geom_20_off + 50 * ccomps * dcomps);

            auto g_xy_0_z_xxy = cbuffer.data(pf_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_z_xxz = cbuffer.data(pf_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_z_xyy = cbuffer.data(pf_geom_20_off + 53 * ccomps * dcomps);

            auto g_xy_0_z_xyz = cbuffer.data(pf_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_z_xzz = cbuffer.data(pf_geom_20_off + 55 * ccomps * dcomps);

            auto g_xy_0_z_yyy = cbuffer.data(pf_geom_20_off + 56 * ccomps * dcomps);

            auto g_xy_0_z_yyz = cbuffer.data(pf_geom_20_off + 57 * ccomps * dcomps);

            auto g_xy_0_z_yzz = cbuffer.data(pf_geom_20_off + 58 * ccomps * dcomps);

            auto g_xy_0_z_zzz = cbuffer.data(pf_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_0_xxx, g_xy_0_0_xxxz, g_xy_0_0_xxy, g_xy_0_0_xxyz, g_xy_0_0_xxz, g_xy_0_0_xxzz, g_xy_0_0_xyy, g_xy_0_0_xyyz, g_xy_0_0_xyz, g_xy_0_0_xyzz, g_xy_0_0_xzz, g_xy_0_0_xzzz, g_xy_0_0_yyy, g_xy_0_0_yyyz, g_xy_0_0_yyz, g_xy_0_0_yyzz, g_xy_0_0_yzz, g_xy_0_0_yzzz, g_xy_0_0_zzz, g_xy_0_0_zzzz, g_xy_0_z_xxx, g_xy_0_z_xxy, g_xy_0_z_xxz, g_xy_0_z_xyy, g_xy_0_z_xyz, g_xy_0_z_xzz, g_xy_0_z_yyy, g_xy_0_z_yyz, g_xy_0_z_yzz, g_xy_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_z_xxx[k] = -g_xy_0_0_xxx[k] * ab_z + g_xy_0_0_xxxz[k];

                g_xy_0_z_xxy[k] = -g_xy_0_0_xxy[k] * ab_z + g_xy_0_0_xxyz[k];

                g_xy_0_z_xxz[k] = -g_xy_0_0_xxz[k] * ab_z + g_xy_0_0_xxzz[k];

                g_xy_0_z_xyy[k] = -g_xy_0_0_xyy[k] * ab_z + g_xy_0_0_xyyz[k];

                g_xy_0_z_xyz[k] = -g_xy_0_0_xyz[k] * ab_z + g_xy_0_0_xyzz[k];

                g_xy_0_z_xzz[k] = -g_xy_0_0_xzz[k] * ab_z + g_xy_0_0_xzzz[k];

                g_xy_0_z_yyy[k] = -g_xy_0_0_yyy[k] * ab_z + g_xy_0_0_yyyz[k];

                g_xy_0_z_yyz[k] = -g_xy_0_0_yyz[k] * ab_z + g_xy_0_0_yyzz[k];

                g_xy_0_z_yzz[k] = -g_xy_0_0_yzz[k] * ab_z + g_xy_0_0_yzzz[k];

                g_xy_0_z_zzz[k] = -g_xy_0_0_zzz[k] * ab_z + g_xy_0_0_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_xz_0_x_xxx = cbuffer.data(pf_geom_20_off + 60 * ccomps * dcomps);

            auto g_xz_0_x_xxy = cbuffer.data(pf_geom_20_off + 61 * ccomps * dcomps);

            auto g_xz_0_x_xxz = cbuffer.data(pf_geom_20_off + 62 * ccomps * dcomps);

            auto g_xz_0_x_xyy = cbuffer.data(pf_geom_20_off + 63 * ccomps * dcomps);

            auto g_xz_0_x_xyz = cbuffer.data(pf_geom_20_off + 64 * ccomps * dcomps);

            auto g_xz_0_x_xzz = cbuffer.data(pf_geom_20_off + 65 * ccomps * dcomps);

            auto g_xz_0_x_yyy = cbuffer.data(pf_geom_20_off + 66 * ccomps * dcomps);

            auto g_xz_0_x_yyz = cbuffer.data(pf_geom_20_off + 67 * ccomps * dcomps);

            auto g_xz_0_x_yzz = cbuffer.data(pf_geom_20_off + 68 * ccomps * dcomps);

            auto g_xz_0_x_zzz = cbuffer.data(pf_geom_20_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_0_xxx, g_xz_0_0_xxxx, g_xz_0_0_xxxy, g_xz_0_0_xxxz, g_xz_0_0_xxy, g_xz_0_0_xxyy, g_xz_0_0_xxyz, g_xz_0_0_xxz, g_xz_0_0_xxzz, g_xz_0_0_xyy, g_xz_0_0_xyyy, g_xz_0_0_xyyz, g_xz_0_0_xyz, g_xz_0_0_xyzz, g_xz_0_0_xzz, g_xz_0_0_xzzz, g_xz_0_0_yyy, g_xz_0_0_yyz, g_xz_0_0_yzz, g_xz_0_0_zzz, g_xz_0_x_xxx, g_xz_0_x_xxy, g_xz_0_x_xxz, g_xz_0_x_xyy, g_xz_0_x_xyz, g_xz_0_x_xzz, g_xz_0_x_yyy, g_xz_0_x_yyz, g_xz_0_x_yzz, g_xz_0_x_zzz, g_z_0_0_xxx, g_z_0_0_xxy, g_z_0_0_xxz, g_z_0_0_xyy, g_z_0_0_xyz, g_z_0_0_xzz, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yzz, g_z_0_0_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_x_xxx[k] = -g_z_0_0_xxx[k] - g_xz_0_0_xxx[k] * ab_x + g_xz_0_0_xxxx[k];

                g_xz_0_x_xxy[k] = -g_z_0_0_xxy[k] - g_xz_0_0_xxy[k] * ab_x + g_xz_0_0_xxxy[k];

                g_xz_0_x_xxz[k] = -g_z_0_0_xxz[k] - g_xz_0_0_xxz[k] * ab_x + g_xz_0_0_xxxz[k];

                g_xz_0_x_xyy[k] = -g_z_0_0_xyy[k] - g_xz_0_0_xyy[k] * ab_x + g_xz_0_0_xxyy[k];

                g_xz_0_x_xyz[k] = -g_z_0_0_xyz[k] - g_xz_0_0_xyz[k] * ab_x + g_xz_0_0_xxyz[k];

                g_xz_0_x_xzz[k] = -g_z_0_0_xzz[k] - g_xz_0_0_xzz[k] * ab_x + g_xz_0_0_xxzz[k];

                g_xz_0_x_yyy[k] = -g_z_0_0_yyy[k] - g_xz_0_0_yyy[k] * ab_x + g_xz_0_0_xyyy[k];

                g_xz_0_x_yyz[k] = -g_z_0_0_yyz[k] - g_xz_0_0_yyz[k] * ab_x + g_xz_0_0_xyyz[k];

                g_xz_0_x_yzz[k] = -g_z_0_0_yzz[k] - g_xz_0_0_yzz[k] * ab_x + g_xz_0_0_xyzz[k];

                g_xz_0_x_zzz[k] = -g_z_0_0_zzz[k] - g_xz_0_0_zzz[k] * ab_x + g_xz_0_0_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_xz_0_y_xxx = cbuffer.data(pf_geom_20_off + 70 * ccomps * dcomps);

            auto g_xz_0_y_xxy = cbuffer.data(pf_geom_20_off + 71 * ccomps * dcomps);

            auto g_xz_0_y_xxz = cbuffer.data(pf_geom_20_off + 72 * ccomps * dcomps);

            auto g_xz_0_y_xyy = cbuffer.data(pf_geom_20_off + 73 * ccomps * dcomps);

            auto g_xz_0_y_xyz = cbuffer.data(pf_geom_20_off + 74 * ccomps * dcomps);

            auto g_xz_0_y_xzz = cbuffer.data(pf_geom_20_off + 75 * ccomps * dcomps);

            auto g_xz_0_y_yyy = cbuffer.data(pf_geom_20_off + 76 * ccomps * dcomps);

            auto g_xz_0_y_yyz = cbuffer.data(pf_geom_20_off + 77 * ccomps * dcomps);

            auto g_xz_0_y_yzz = cbuffer.data(pf_geom_20_off + 78 * ccomps * dcomps);

            auto g_xz_0_y_zzz = cbuffer.data(pf_geom_20_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_0_xxx, g_xz_0_0_xxxy, g_xz_0_0_xxy, g_xz_0_0_xxyy, g_xz_0_0_xxyz, g_xz_0_0_xxz, g_xz_0_0_xyy, g_xz_0_0_xyyy, g_xz_0_0_xyyz, g_xz_0_0_xyz, g_xz_0_0_xyzz, g_xz_0_0_xzz, g_xz_0_0_yyy, g_xz_0_0_yyyy, g_xz_0_0_yyyz, g_xz_0_0_yyz, g_xz_0_0_yyzz, g_xz_0_0_yzz, g_xz_0_0_yzzz, g_xz_0_0_zzz, g_xz_0_y_xxx, g_xz_0_y_xxy, g_xz_0_y_xxz, g_xz_0_y_xyy, g_xz_0_y_xyz, g_xz_0_y_xzz, g_xz_0_y_yyy, g_xz_0_y_yyz, g_xz_0_y_yzz, g_xz_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_y_xxx[k] = -g_xz_0_0_xxx[k] * ab_y + g_xz_0_0_xxxy[k];

                g_xz_0_y_xxy[k] = -g_xz_0_0_xxy[k] * ab_y + g_xz_0_0_xxyy[k];

                g_xz_0_y_xxz[k] = -g_xz_0_0_xxz[k] * ab_y + g_xz_0_0_xxyz[k];

                g_xz_0_y_xyy[k] = -g_xz_0_0_xyy[k] * ab_y + g_xz_0_0_xyyy[k];

                g_xz_0_y_xyz[k] = -g_xz_0_0_xyz[k] * ab_y + g_xz_0_0_xyyz[k];

                g_xz_0_y_xzz[k] = -g_xz_0_0_xzz[k] * ab_y + g_xz_0_0_xyzz[k];

                g_xz_0_y_yyy[k] = -g_xz_0_0_yyy[k] * ab_y + g_xz_0_0_yyyy[k];

                g_xz_0_y_yyz[k] = -g_xz_0_0_yyz[k] * ab_y + g_xz_0_0_yyyz[k];

                g_xz_0_y_yzz[k] = -g_xz_0_0_yzz[k] * ab_y + g_xz_0_0_yyzz[k];

                g_xz_0_y_zzz[k] = -g_xz_0_0_zzz[k] * ab_y + g_xz_0_0_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_xz_0_z_xxx = cbuffer.data(pf_geom_20_off + 80 * ccomps * dcomps);

            auto g_xz_0_z_xxy = cbuffer.data(pf_geom_20_off + 81 * ccomps * dcomps);

            auto g_xz_0_z_xxz = cbuffer.data(pf_geom_20_off + 82 * ccomps * dcomps);

            auto g_xz_0_z_xyy = cbuffer.data(pf_geom_20_off + 83 * ccomps * dcomps);

            auto g_xz_0_z_xyz = cbuffer.data(pf_geom_20_off + 84 * ccomps * dcomps);

            auto g_xz_0_z_xzz = cbuffer.data(pf_geom_20_off + 85 * ccomps * dcomps);

            auto g_xz_0_z_yyy = cbuffer.data(pf_geom_20_off + 86 * ccomps * dcomps);

            auto g_xz_0_z_yyz = cbuffer.data(pf_geom_20_off + 87 * ccomps * dcomps);

            auto g_xz_0_z_yzz = cbuffer.data(pf_geom_20_off + 88 * ccomps * dcomps);

            auto g_xz_0_z_zzz = cbuffer.data(pf_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxx, g_x_0_0_xxy, g_x_0_0_xxz, g_x_0_0_xyy, g_x_0_0_xyz, g_x_0_0_xzz, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yzz, g_x_0_0_zzz, g_xz_0_0_xxx, g_xz_0_0_xxxz, g_xz_0_0_xxy, g_xz_0_0_xxyz, g_xz_0_0_xxz, g_xz_0_0_xxzz, g_xz_0_0_xyy, g_xz_0_0_xyyz, g_xz_0_0_xyz, g_xz_0_0_xyzz, g_xz_0_0_xzz, g_xz_0_0_xzzz, g_xz_0_0_yyy, g_xz_0_0_yyyz, g_xz_0_0_yyz, g_xz_0_0_yyzz, g_xz_0_0_yzz, g_xz_0_0_yzzz, g_xz_0_0_zzz, g_xz_0_0_zzzz, g_xz_0_z_xxx, g_xz_0_z_xxy, g_xz_0_z_xxz, g_xz_0_z_xyy, g_xz_0_z_xyz, g_xz_0_z_xzz, g_xz_0_z_yyy, g_xz_0_z_yyz, g_xz_0_z_yzz, g_xz_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_z_xxx[k] = -g_x_0_0_xxx[k] - g_xz_0_0_xxx[k] * ab_z + g_xz_0_0_xxxz[k];

                g_xz_0_z_xxy[k] = -g_x_0_0_xxy[k] - g_xz_0_0_xxy[k] * ab_z + g_xz_0_0_xxyz[k];

                g_xz_0_z_xxz[k] = -g_x_0_0_xxz[k] - g_xz_0_0_xxz[k] * ab_z + g_xz_0_0_xxzz[k];

                g_xz_0_z_xyy[k] = -g_x_0_0_xyy[k] - g_xz_0_0_xyy[k] * ab_z + g_xz_0_0_xyyz[k];

                g_xz_0_z_xyz[k] = -g_x_0_0_xyz[k] - g_xz_0_0_xyz[k] * ab_z + g_xz_0_0_xyzz[k];

                g_xz_0_z_xzz[k] = -g_x_0_0_xzz[k] - g_xz_0_0_xzz[k] * ab_z + g_xz_0_0_xzzz[k];

                g_xz_0_z_yyy[k] = -g_x_0_0_yyy[k] - g_xz_0_0_yyy[k] * ab_z + g_xz_0_0_yyyz[k];

                g_xz_0_z_yyz[k] = -g_x_0_0_yyz[k] - g_xz_0_0_yyz[k] * ab_z + g_xz_0_0_yyzz[k];

                g_xz_0_z_yzz[k] = -g_x_0_0_yzz[k] - g_xz_0_0_yzz[k] * ab_z + g_xz_0_0_yzzz[k];

                g_xz_0_z_zzz[k] = -g_x_0_0_zzz[k] - g_xz_0_0_zzz[k] * ab_z + g_xz_0_0_zzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_yy_0_x_xxx = cbuffer.data(pf_geom_20_off + 90 * ccomps * dcomps);

            auto g_yy_0_x_xxy = cbuffer.data(pf_geom_20_off + 91 * ccomps * dcomps);

            auto g_yy_0_x_xxz = cbuffer.data(pf_geom_20_off + 92 * ccomps * dcomps);

            auto g_yy_0_x_xyy = cbuffer.data(pf_geom_20_off + 93 * ccomps * dcomps);

            auto g_yy_0_x_xyz = cbuffer.data(pf_geom_20_off + 94 * ccomps * dcomps);

            auto g_yy_0_x_xzz = cbuffer.data(pf_geom_20_off + 95 * ccomps * dcomps);

            auto g_yy_0_x_yyy = cbuffer.data(pf_geom_20_off + 96 * ccomps * dcomps);

            auto g_yy_0_x_yyz = cbuffer.data(pf_geom_20_off + 97 * ccomps * dcomps);

            auto g_yy_0_x_yzz = cbuffer.data(pf_geom_20_off + 98 * ccomps * dcomps);

            auto g_yy_0_x_zzz = cbuffer.data(pf_geom_20_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_0_xxx, g_yy_0_0_xxxx, g_yy_0_0_xxxy, g_yy_0_0_xxxz, g_yy_0_0_xxy, g_yy_0_0_xxyy, g_yy_0_0_xxyz, g_yy_0_0_xxz, g_yy_0_0_xxzz, g_yy_0_0_xyy, g_yy_0_0_xyyy, g_yy_0_0_xyyz, g_yy_0_0_xyz, g_yy_0_0_xyzz, g_yy_0_0_xzz, g_yy_0_0_xzzz, g_yy_0_0_yyy, g_yy_0_0_yyz, g_yy_0_0_yzz, g_yy_0_0_zzz, g_yy_0_x_xxx, g_yy_0_x_xxy, g_yy_0_x_xxz, g_yy_0_x_xyy, g_yy_0_x_xyz, g_yy_0_x_xzz, g_yy_0_x_yyy, g_yy_0_x_yyz, g_yy_0_x_yzz, g_yy_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_x_xxx[k] = -g_yy_0_0_xxx[k] * ab_x + g_yy_0_0_xxxx[k];

                g_yy_0_x_xxy[k] = -g_yy_0_0_xxy[k] * ab_x + g_yy_0_0_xxxy[k];

                g_yy_0_x_xxz[k] = -g_yy_0_0_xxz[k] * ab_x + g_yy_0_0_xxxz[k];

                g_yy_0_x_xyy[k] = -g_yy_0_0_xyy[k] * ab_x + g_yy_0_0_xxyy[k];

                g_yy_0_x_xyz[k] = -g_yy_0_0_xyz[k] * ab_x + g_yy_0_0_xxyz[k];

                g_yy_0_x_xzz[k] = -g_yy_0_0_xzz[k] * ab_x + g_yy_0_0_xxzz[k];

                g_yy_0_x_yyy[k] = -g_yy_0_0_yyy[k] * ab_x + g_yy_0_0_xyyy[k];

                g_yy_0_x_yyz[k] = -g_yy_0_0_yyz[k] * ab_x + g_yy_0_0_xyyz[k];

                g_yy_0_x_yzz[k] = -g_yy_0_0_yzz[k] * ab_x + g_yy_0_0_xyzz[k];

                g_yy_0_x_zzz[k] = -g_yy_0_0_zzz[k] * ab_x + g_yy_0_0_xzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_yy_0_y_xxx = cbuffer.data(pf_geom_20_off + 100 * ccomps * dcomps);

            auto g_yy_0_y_xxy = cbuffer.data(pf_geom_20_off + 101 * ccomps * dcomps);

            auto g_yy_0_y_xxz = cbuffer.data(pf_geom_20_off + 102 * ccomps * dcomps);

            auto g_yy_0_y_xyy = cbuffer.data(pf_geom_20_off + 103 * ccomps * dcomps);

            auto g_yy_0_y_xyz = cbuffer.data(pf_geom_20_off + 104 * ccomps * dcomps);

            auto g_yy_0_y_xzz = cbuffer.data(pf_geom_20_off + 105 * ccomps * dcomps);

            auto g_yy_0_y_yyy = cbuffer.data(pf_geom_20_off + 106 * ccomps * dcomps);

            auto g_yy_0_y_yyz = cbuffer.data(pf_geom_20_off + 107 * ccomps * dcomps);

            auto g_yy_0_y_yzz = cbuffer.data(pf_geom_20_off + 108 * ccomps * dcomps);

            auto g_yy_0_y_zzz = cbuffer.data(pf_geom_20_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxx, g_y_0_0_xxy, g_y_0_0_xxz, g_y_0_0_xyy, g_y_0_0_xyz, g_y_0_0_xzz, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yzz, g_y_0_0_zzz, g_yy_0_0_xxx, g_yy_0_0_xxxy, g_yy_0_0_xxy, g_yy_0_0_xxyy, g_yy_0_0_xxyz, g_yy_0_0_xxz, g_yy_0_0_xyy, g_yy_0_0_xyyy, g_yy_0_0_xyyz, g_yy_0_0_xyz, g_yy_0_0_xyzz, g_yy_0_0_xzz, g_yy_0_0_yyy, g_yy_0_0_yyyy, g_yy_0_0_yyyz, g_yy_0_0_yyz, g_yy_0_0_yyzz, g_yy_0_0_yzz, g_yy_0_0_yzzz, g_yy_0_0_zzz, g_yy_0_y_xxx, g_yy_0_y_xxy, g_yy_0_y_xxz, g_yy_0_y_xyy, g_yy_0_y_xyz, g_yy_0_y_xzz, g_yy_0_y_yyy, g_yy_0_y_yyz, g_yy_0_y_yzz, g_yy_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_y_xxx[k] = -2.0 * g_y_0_0_xxx[k] - g_yy_0_0_xxx[k] * ab_y + g_yy_0_0_xxxy[k];

                g_yy_0_y_xxy[k] = -2.0 * g_y_0_0_xxy[k] - g_yy_0_0_xxy[k] * ab_y + g_yy_0_0_xxyy[k];

                g_yy_0_y_xxz[k] = -2.0 * g_y_0_0_xxz[k] - g_yy_0_0_xxz[k] * ab_y + g_yy_0_0_xxyz[k];

                g_yy_0_y_xyy[k] = -2.0 * g_y_0_0_xyy[k] - g_yy_0_0_xyy[k] * ab_y + g_yy_0_0_xyyy[k];

                g_yy_0_y_xyz[k] = -2.0 * g_y_0_0_xyz[k] - g_yy_0_0_xyz[k] * ab_y + g_yy_0_0_xyyz[k];

                g_yy_0_y_xzz[k] = -2.0 * g_y_0_0_xzz[k] - g_yy_0_0_xzz[k] * ab_y + g_yy_0_0_xyzz[k];

                g_yy_0_y_yyy[k] = -2.0 * g_y_0_0_yyy[k] - g_yy_0_0_yyy[k] * ab_y + g_yy_0_0_yyyy[k];

                g_yy_0_y_yyz[k] = -2.0 * g_y_0_0_yyz[k] - g_yy_0_0_yyz[k] * ab_y + g_yy_0_0_yyyz[k];

                g_yy_0_y_yzz[k] = -2.0 * g_y_0_0_yzz[k] - g_yy_0_0_yzz[k] * ab_y + g_yy_0_0_yyzz[k];

                g_yy_0_y_zzz[k] = -2.0 * g_y_0_0_zzz[k] - g_yy_0_0_zzz[k] * ab_y + g_yy_0_0_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_yy_0_z_xxx = cbuffer.data(pf_geom_20_off + 110 * ccomps * dcomps);

            auto g_yy_0_z_xxy = cbuffer.data(pf_geom_20_off + 111 * ccomps * dcomps);

            auto g_yy_0_z_xxz = cbuffer.data(pf_geom_20_off + 112 * ccomps * dcomps);

            auto g_yy_0_z_xyy = cbuffer.data(pf_geom_20_off + 113 * ccomps * dcomps);

            auto g_yy_0_z_xyz = cbuffer.data(pf_geom_20_off + 114 * ccomps * dcomps);

            auto g_yy_0_z_xzz = cbuffer.data(pf_geom_20_off + 115 * ccomps * dcomps);

            auto g_yy_0_z_yyy = cbuffer.data(pf_geom_20_off + 116 * ccomps * dcomps);

            auto g_yy_0_z_yyz = cbuffer.data(pf_geom_20_off + 117 * ccomps * dcomps);

            auto g_yy_0_z_yzz = cbuffer.data(pf_geom_20_off + 118 * ccomps * dcomps);

            auto g_yy_0_z_zzz = cbuffer.data(pf_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_0_xxx, g_yy_0_0_xxxz, g_yy_0_0_xxy, g_yy_0_0_xxyz, g_yy_0_0_xxz, g_yy_0_0_xxzz, g_yy_0_0_xyy, g_yy_0_0_xyyz, g_yy_0_0_xyz, g_yy_0_0_xyzz, g_yy_0_0_xzz, g_yy_0_0_xzzz, g_yy_0_0_yyy, g_yy_0_0_yyyz, g_yy_0_0_yyz, g_yy_0_0_yyzz, g_yy_0_0_yzz, g_yy_0_0_yzzz, g_yy_0_0_zzz, g_yy_0_0_zzzz, g_yy_0_z_xxx, g_yy_0_z_xxy, g_yy_0_z_xxz, g_yy_0_z_xyy, g_yy_0_z_xyz, g_yy_0_z_xzz, g_yy_0_z_yyy, g_yy_0_z_yyz, g_yy_0_z_yzz, g_yy_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_z_xxx[k] = -g_yy_0_0_xxx[k] * ab_z + g_yy_0_0_xxxz[k];

                g_yy_0_z_xxy[k] = -g_yy_0_0_xxy[k] * ab_z + g_yy_0_0_xxyz[k];

                g_yy_0_z_xxz[k] = -g_yy_0_0_xxz[k] * ab_z + g_yy_0_0_xxzz[k];

                g_yy_0_z_xyy[k] = -g_yy_0_0_xyy[k] * ab_z + g_yy_0_0_xyyz[k];

                g_yy_0_z_xyz[k] = -g_yy_0_0_xyz[k] * ab_z + g_yy_0_0_xyzz[k];

                g_yy_0_z_xzz[k] = -g_yy_0_0_xzz[k] * ab_z + g_yy_0_0_xzzz[k];

                g_yy_0_z_yyy[k] = -g_yy_0_0_yyy[k] * ab_z + g_yy_0_0_yyyz[k];

                g_yy_0_z_yyz[k] = -g_yy_0_0_yyz[k] * ab_z + g_yy_0_0_yyzz[k];

                g_yy_0_z_yzz[k] = -g_yy_0_0_yzz[k] * ab_z + g_yy_0_0_yzzz[k];

                g_yy_0_z_zzz[k] = -g_yy_0_0_zzz[k] * ab_z + g_yy_0_0_zzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_yz_0_x_xxx = cbuffer.data(pf_geom_20_off + 120 * ccomps * dcomps);

            auto g_yz_0_x_xxy = cbuffer.data(pf_geom_20_off + 121 * ccomps * dcomps);

            auto g_yz_0_x_xxz = cbuffer.data(pf_geom_20_off + 122 * ccomps * dcomps);

            auto g_yz_0_x_xyy = cbuffer.data(pf_geom_20_off + 123 * ccomps * dcomps);

            auto g_yz_0_x_xyz = cbuffer.data(pf_geom_20_off + 124 * ccomps * dcomps);

            auto g_yz_0_x_xzz = cbuffer.data(pf_geom_20_off + 125 * ccomps * dcomps);

            auto g_yz_0_x_yyy = cbuffer.data(pf_geom_20_off + 126 * ccomps * dcomps);

            auto g_yz_0_x_yyz = cbuffer.data(pf_geom_20_off + 127 * ccomps * dcomps);

            auto g_yz_0_x_yzz = cbuffer.data(pf_geom_20_off + 128 * ccomps * dcomps);

            auto g_yz_0_x_zzz = cbuffer.data(pf_geom_20_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_0_xxx, g_yz_0_0_xxxx, g_yz_0_0_xxxy, g_yz_0_0_xxxz, g_yz_0_0_xxy, g_yz_0_0_xxyy, g_yz_0_0_xxyz, g_yz_0_0_xxz, g_yz_0_0_xxzz, g_yz_0_0_xyy, g_yz_0_0_xyyy, g_yz_0_0_xyyz, g_yz_0_0_xyz, g_yz_0_0_xyzz, g_yz_0_0_xzz, g_yz_0_0_xzzz, g_yz_0_0_yyy, g_yz_0_0_yyz, g_yz_0_0_yzz, g_yz_0_0_zzz, g_yz_0_x_xxx, g_yz_0_x_xxy, g_yz_0_x_xxz, g_yz_0_x_xyy, g_yz_0_x_xyz, g_yz_0_x_xzz, g_yz_0_x_yyy, g_yz_0_x_yyz, g_yz_0_x_yzz, g_yz_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_x_xxx[k] = -g_yz_0_0_xxx[k] * ab_x + g_yz_0_0_xxxx[k];

                g_yz_0_x_xxy[k] = -g_yz_0_0_xxy[k] * ab_x + g_yz_0_0_xxxy[k];

                g_yz_0_x_xxz[k] = -g_yz_0_0_xxz[k] * ab_x + g_yz_0_0_xxxz[k];

                g_yz_0_x_xyy[k] = -g_yz_0_0_xyy[k] * ab_x + g_yz_0_0_xxyy[k];

                g_yz_0_x_xyz[k] = -g_yz_0_0_xyz[k] * ab_x + g_yz_0_0_xxyz[k];

                g_yz_0_x_xzz[k] = -g_yz_0_0_xzz[k] * ab_x + g_yz_0_0_xxzz[k];

                g_yz_0_x_yyy[k] = -g_yz_0_0_yyy[k] * ab_x + g_yz_0_0_xyyy[k];

                g_yz_0_x_yyz[k] = -g_yz_0_0_yyz[k] * ab_x + g_yz_0_0_xyyz[k];

                g_yz_0_x_yzz[k] = -g_yz_0_0_yzz[k] * ab_x + g_yz_0_0_xyzz[k];

                g_yz_0_x_zzz[k] = -g_yz_0_0_zzz[k] * ab_x + g_yz_0_0_xzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_yz_0_y_xxx = cbuffer.data(pf_geom_20_off + 130 * ccomps * dcomps);

            auto g_yz_0_y_xxy = cbuffer.data(pf_geom_20_off + 131 * ccomps * dcomps);

            auto g_yz_0_y_xxz = cbuffer.data(pf_geom_20_off + 132 * ccomps * dcomps);

            auto g_yz_0_y_xyy = cbuffer.data(pf_geom_20_off + 133 * ccomps * dcomps);

            auto g_yz_0_y_xyz = cbuffer.data(pf_geom_20_off + 134 * ccomps * dcomps);

            auto g_yz_0_y_xzz = cbuffer.data(pf_geom_20_off + 135 * ccomps * dcomps);

            auto g_yz_0_y_yyy = cbuffer.data(pf_geom_20_off + 136 * ccomps * dcomps);

            auto g_yz_0_y_yyz = cbuffer.data(pf_geom_20_off + 137 * ccomps * dcomps);

            auto g_yz_0_y_yzz = cbuffer.data(pf_geom_20_off + 138 * ccomps * dcomps);

            auto g_yz_0_y_zzz = cbuffer.data(pf_geom_20_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_0_xxx, g_yz_0_0_xxxy, g_yz_0_0_xxy, g_yz_0_0_xxyy, g_yz_0_0_xxyz, g_yz_0_0_xxz, g_yz_0_0_xyy, g_yz_0_0_xyyy, g_yz_0_0_xyyz, g_yz_0_0_xyz, g_yz_0_0_xyzz, g_yz_0_0_xzz, g_yz_0_0_yyy, g_yz_0_0_yyyy, g_yz_0_0_yyyz, g_yz_0_0_yyz, g_yz_0_0_yyzz, g_yz_0_0_yzz, g_yz_0_0_yzzz, g_yz_0_0_zzz, g_yz_0_y_xxx, g_yz_0_y_xxy, g_yz_0_y_xxz, g_yz_0_y_xyy, g_yz_0_y_xyz, g_yz_0_y_xzz, g_yz_0_y_yyy, g_yz_0_y_yyz, g_yz_0_y_yzz, g_yz_0_y_zzz, g_z_0_0_xxx, g_z_0_0_xxy, g_z_0_0_xxz, g_z_0_0_xyy, g_z_0_0_xyz, g_z_0_0_xzz, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yzz, g_z_0_0_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_y_xxx[k] = -g_z_0_0_xxx[k] - g_yz_0_0_xxx[k] * ab_y + g_yz_0_0_xxxy[k];

                g_yz_0_y_xxy[k] = -g_z_0_0_xxy[k] - g_yz_0_0_xxy[k] * ab_y + g_yz_0_0_xxyy[k];

                g_yz_0_y_xxz[k] = -g_z_0_0_xxz[k] - g_yz_0_0_xxz[k] * ab_y + g_yz_0_0_xxyz[k];

                g_yz_0_y_xyy[k] = -g_z_0_0_xyy[k] - g_yz_0_0_xyy[k] * ab_y + g_yz_0_0_xyyy[k];

                g_yz_0_y_xyz[k] = -g_z_0_0_xyz[k] - g_yz_0_0_xyz[k] * ab_y + g_yz_0_0_xyyz[k];

                g_yz_0_y_xzz[k] = -g_z_0_0_xzz[k] - g_yz_0_0_xzz[k] * ab_y + g_yz_0_0_xyzz[k];

                g_yz_0_y_yyy[k] = -g_z_0_0_yyy[k] - g_yz_0_0_yyy[k] * ab_y + g_yz_0_0_yyyy[k];

                g_yz_0_y_yyz[k] = -g_z_0_0_yyz[k] - g_yz_0_0_yyz[k] * ab_y + g_yz_0_0_yyyz[k];

                g_yz_0_y_yzz[k] = -g_z_0_0_yzz[k] - g_yz_0_0_yzz[k] * ab_y + g_yz_0_0_yyzz[k];

                g_yz_0_y_zzz[k] = -g_z_0_0_zzz[k] - g_yz_0_0_zzz[k] * ab_y + g_yz_0_0_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_yz_0_z_xxx = cbuffer.data(pf_geom_20_off + 140 * ccomps * dcomps);

            auto g_yz_0_z_xxy = cbuffer.data(pf_geom_20_off + 141 * ccomps * dcomps);

            auto g_yz_0_z_xxz = cbuffer.data(pf_geom_20_off + 142 * ccomps * dcomps);

            auto g_yz_0_z_xyy = cbuffer.data(pf_geom_20_off + 143 * ccomps * dcomps);

            auto g_yz_0_z_xyz = cbuffer.data(pf_geom_20_off + 144 * ccomps * dcomps);

            auto g_yz_0_z_xzz = cbuffer.data(pf_geom_20_off + 145 * ccomps * dcomps);

            auto g_yz_0_z_yyy = cbuffer.data(pf_geom_20_off + 146 * ccomps * dcomps);

            auto g_yz_0_z_yyz = cbuffer.data(pf_geom_20_off + 147 * ccomps * dcomps);

            auto g_yz_0_z_yzz = cbuffer.data(pf_geom_20_off + 148 * ccomps * dcomps);

            auto g_yz_0_z_zzz = cbuffer.data(pf_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxx, g_y_0_0_xxy, g_y_0_0_xxz, g_y_0_0_xyy, g_y_0_0_xyz, g_y_0_0_xzz, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yzz, g_y_0_0_zzz, g_yz_0_0_xxx, g_yz_0_0_xxxz, g_yz_0_0_xxy, g_yz_0_0_xxyz, g_yz_0_0_xxz, g_yz_0_0_xxzz, g_yz_0_0_xyy, g_yz_0_0_xyyz, g_yz_0_0_xyz, g_yz_0_0_xyzz, g_yz_0_0_xzz, g_yz_0_0_xzzz, g_yz_0_0_yyy, g_yz_0_0_yyyz, g_yz_0_0_yyz, g_yz_0_0_yyzz, g_yz_0_0_yzz, g_yz_0_0_yzzz, g_yz_0_0_zzz, g_yz_0_0_zzzz, g_yz_0_z_xxx, g_yz_0_z_xxy, g_yz_0_z_xxz, g_yz_0_z_xyy, g_yz_0_z_xyz, g_yz_0_z_xzz, g_yz_0_z_yyy, g_yz_0_z_yyz, g_yz_0_z_yzz, g_yz_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_z_xxx[k] = -g_y_0_0_xxx[k] - g_yz_0_0_xxx[k] * ab_z + g_yz_0_0_xxxz[k];

                g_yz_0_z_xxy[k] = -g_y_0_0_xxy[k] - g_yz_0_0_xxy[k] * ab_z + g_yz_0_0_xxyz[k];

                g_yz_0_z_xxz[k] = -g_y_0_0_xxz[k] - g_yz_0_0_xxz[k] * ab_z + g_yz_0_0_xxzz[k];

                g_yz_0_z_xyy[k] = -g_y_0_0_xyy[k] - g_yz_0_0_xyy[k] * ab_z + g_yz_0_0_xyyz[k];

                g_yz_0_z_xyz[k] = -g_y_0_0_xyz[k] - g_yz_0_0_xyz[k] * ab_z + g_yz_0_0_xyzz[k];

                g_yz_0_z_xzz[k] = -g_y_0_0_xzz[k] - g_yz_0_0_xzz[k] * ab_z + g_yz_0_0_xzzz[k];

                g_yz_0_z_yyy[k] = -g_y_0_0_yyy[k] - g_yz_0_0_yyy[k] * ab_z + g_yz_0_0_yyyz[k];

                g_yz_0_z_yyz[k] = -g_y_0_0_yyz[k] - g_yz_0_0_yyz[k] * ab_z + g_yz_0_0_yyzz[k];

                g_yz_0_z_yzz[k] = -g_y_0_0_yzz[k] - g_yz_0_0_yzz[k] * ab_z + g_yz_0_0_yzzz[k];

                g_yz_0_z_zzz[k] = -g_y_0_0_zzz[k] - g_yz_0_0_zzz[k] * ab_z + g_yz_0_0_zzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_zz_0_x_xxx = cbuffer.data(pf_geom_20_off + 150 * ccomps * dcomps);

            auto g_zz_0_x_xxy = cbuffer.data(pf_geom_20_off + 151 * ccomps * dcomps);

            auto g_zz_0_x_xxz = cbuffer.data(pf_geom_20_off + 152 * ccomps * dcomps);

            auto g_zz_0_x_xyy = cbuffer.data(pf_geom_20_off + 153 * ccomps * dcomps);

            auto g_zz_0_x_xyz = cbuffer.data(pf_geom_20_off + 154 * ccomps * dcomps);

            auto g_zz_0_x_xzz = cbuffer.data(pf_geom_20_off + 155 * ccomps * dcomps);

            auto g_zz_0_x_yyy = cbuffer.data(pf_geom_20_off + 156 * ccomps * dcomps);

            auto g_zz_0_x_yyz = cbuffer.data(pf_geom_20_off + 157 * ccomps * dcomps);

            auto g_zz_0_x_yzz = cbuffer.data(pf_geom_20_off + 158 * ccomps * dcomps);

            auto g_zz_0_x_zzz = cbuffer.data(pf_geom_20_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_0_xxx, g_zz_0_0_xxxx, g_zz_0_0_xxxy, g_zz_0_0_xxxz, g_zz_0_0_xxy, g_zz_0_0_xxyy, g_zz_0_0_xxyz, g_zz_0_0_xxz, g_zz_0_0_xxzz, g_zz_0_0_xyy, g_zz_0_0_xyyy, g_zz_0_0_xyyz, g_zz_0_0_xyz, g_zz_0_0_xyzz, g_zz_0_0_xzz, g_zz_0_0_xzzz, g_zz_0_0_yyy, g_zz_0_0_yyz, g_zz_0_0_yzz, g_zz_0_0_zzz, g_zz_0_x_xxx, g_zz_0_x_xxy, g_zz_0_x_xxz, g_zz_0_x_xyy, g_zz_0_x_xyz, g_zz_0_x_xzz, g_zz_0_x_yyy, g_zz_0_x_yyz, g_zz_0_x_yzz, g_zz_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_x_xxx[k] = -g_zz_0_0_xxx[k] * ab_x + g_zz_0_0_xxxx[k];

                g_zz_0_x_xxy[k] = -g_zz_0_0_xxy[k] * ab_x + g_zz_0_0_xxxy[k];

                g_zz_0_x_xxz[k] = -g_zz_0_0_xxz[k] * ab_x + g_zz_0_0_xxxz[k];

                g_zz_0_x_xyy[k] = -g_zz_0_0_xyy[k] * ab_x + g_zz_0_0_xxyy[k];

                g_zz_0_x_xyz[k] = -g_zz_0_0_xyz[k] * ab_x + g_zz_0_0_xxyz[k];

                g_zz_0_x_xzz[k] = -g_zz_0_0_xzz[k] * ab_x + g_zz_0_0_xxzz[k];

                g_zz_0_x_yyy[k] = -g_zz_0_0_yyy[k] * ab_x + g_zz_0_0_xyyy[k];

                g_zz_0_x_yyz[k] = -g_zz_0_0_yyz[k] * ab_x + g_zz_0_0_xyyz[k];

                g_zz_0_x_yzz[k] = -g_zz_0_0_yzz[k] * ab_x + g_zz_0_0_xyzz[k];

                g_zz_0_x_zzz[k] = -g_zz_0_0_zzz[k] * ab_x + g_zz_0_0_xzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_zz_0_y_xxx = cbuffer.data(pf_geom_20_off + 160 * ccomps * dcomps);

            auto g_zz_0_y_xxy = cbuffer.data(pf_geom_20_off + 161 * ccomps * dcomps);

            auto g_zz_0_y_xxz = cbuffer.data(pf_geom_20_off + 162 * ccomps * dcomps);

            auto g_zz_0_y_xyy = cbuffer.data(pf_geom_20_off + 163 * ccomps * dcomps);

            auto g_zz_0_y_xyz = cbuffer.data(pf_geom_20_off + 164 * ccomps * dcomps);

            auto g_zz_0_y_xzz = cbuffer.data(pf_geom_20_off + 165 * ccomps * dcomps);

            auto g_zz_0_y_yyy = cbuffer.data(pf_geom_20_off + 166 * ccomps * dcomps);

            auto g_zz_0_y_yyz = cbuffer.data(pf_geom_20_off + 167 * ccomps * dcomps);

            auto g_zz_0_y_yzz = cbuffer.data(pf_geom_20_off + 168 * ccomps * dcomps);

            auto g_zz_0_y_zzz = cbuffer.data(pf_geom_20_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_0_xxx, g_zz_0_0_xxxy, g_zz_0_0_xxy, g_zz_0_0_xxyy, g_zz_0_0_xxyz, g_zz_0_0_xxz, g_zz_0_0_xyy, g_zz_0_0_xyyy, g_zz_0_0_xyyz, g_zz_0_0_xyz, g_zz_0_0_xyzz, g_zz_0_0_xzz, g_zz_0_0_yyy, g_zz_0_0_yyyy, g_zz_0_0_yyyz, g_zz_0_0_yyz, g_zz_0_0_yyzz, g_zz_0_0_yzz, g_zz_0_0_yzzz, g_zz_0_0_zzz, g_zz_0_y_xxx, g_zz_0_y_xxy, g_zz_0_y_xxz, g_zz_0_y_xyy, g_zz_0_y_xyz, g_zz_0_y_xzz, g_zz_0_y_yyy, g_zz_0_y_yyz, g_zz_0_y_yzz, g_zz_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_y_xxx[k] = -g_zz_0_0_xxx[k] * ab_y + g_zz_0_0_xxxy[k];

                g_zz_0_y_xxy[k] = -g_zz_0_0_xxy[k] * ab_y + g_zz_0_0_xxyy[k];

                g_zz_0_y_xxz[k] = -g_zz_0_0_xxz[k] * ab_y + g_zz_0_0_xxyz[k];

                g_zz_0_y_xyy[k] = -g_zz_0_0_xyy[k] * ab_y + g_zz_0_0_xyyy[k];

                g_zz_0_y_xyz[k] = -g_zz_0_0_xyz[k] * ab_y + g_zz_0_0_xyyz[k];

                g_zz_0_y_xzz[k] = -g_zz_0_0_xzz[k] * ab_y + g_zz_0_0_xyzz[k];

                g_zz_0_y_yyy[k] = -g_zz_0_0_yyy[k] * ab_y + g_zz_0_0_yyyy[k];

                g_zz_0_y_yyz[k] = -g_zz_0_0_yyz[k] * ab_y + g_zz_0_0_yyyz[k];

                g_zz_0_y_yzz[k] = -g_zz_0_0_yzz[k] * ab_y + g_zz_0_0_yyzz[k];

                g_zz_0_y_zzz[k] = -g_zz_0_0_zzz[k] * ab_y + g_zz_0_0_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_zz_0_z_xxx = cbuffer.data(pf_geom_20_off + 170 * ccomps * dcomps);

            auto g_zz_0_z_xxy = cbuffer.data(pf_geom_20_off + 171 * ccomps * dcomps);

            auto g_zz_0_z_xxz = cbuffer.data(pf_geom_20_off + 172 * ccomps * dcomps);

            auto g_zz_0_z_xyy = cbuffer.data(pf_geom_20_off + 173 * ccomps * dcomps);

            auto g_zz_0_z_xyz = cbuffer.data(pf_geom_20_off + 174 * ccomps * dcomps);

            auto g_zz_0_z_xzz = cbuffer.data(pf_geom_20_off + 175 * ccomps * dcomps);

            auto g_zz_0_z_yyy = cbuffer.data(pf_geom_20_off + 176 * ccomps * dcomps);

            auto g_zz_0_z_yyz = cbuffer.data(pf_geom_20_off + 177 * ccomps * dcomps);

            auto g_zz_0_z_yzz = cbuffer.data(pf_geom_20_off + 178 * ccomps * dcomps);

            auto g_zz_0_z_zzz = cbuffer.data(pf_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxx, g_z_0_0_xxy, g_z_0_0_xxz, g_z_0_0_xyy, g_z_0_0_xyz, g_z_0_0_xzz, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yzz, g_z_0_0_zzz, g_zz_0_0_xxx, g_zz_0_0_xxxz, g_zz_0_0_xxy, g_zz_0_0_xxyz, g_zz_0_0_xxz, g_zz_0_0_xxzz, g_zz_0_0_xyy, g_zz_0_0_xyyz, g_zz_0_0_xyz, g_zz_0_0_xyzz, g_zz_0_0_xzz, g_zz_0_0_xzzz, g_zz_0_0_yyy, g_zz_0_0_yyyz, g_zz_0_0_yyz, g_zz_0_0_yyzz, g_zz_0_0_yzz, g_zz_0_0_yzzz, g_zz_0_0_zzz, g_zz_0_0_zzzz, g_zz_0_z_xxx, g_zz_0_z_xxy, g_zz_0_z_xxz, g_zz_0_z_xyy, g_zz_0_z_xyz, g_zz_0_z_xzz, g_zz_0_z_yyy, g_zz_0_z_yyz, g_zz_0_z_yzz, g_zz_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_z_xxx[k] = -2.0 * g_z_0_0_xxx[k] - g_zz_0_0_xxx[k] * ab_z + g_zz_0_0_xxxz[k];

                g_zz_0_z_xxy[k] = -2.0 * g_z_0_0_xxy[k] - g_zz_0_0_xxy[k] * ab_z + g_zz_0_0_xxyz[k];

                g_zz_0_z_xxz[k] = -2.0 * g_z_0_0_xxz[k] - g_zz_0_0_xxz[k] * ab_z + g_zz_0_0_xxzz[k];

                g_zz_0_z_xyy[k] = -2.0 * g_z_0_0_xyy[k] - g_zz_0_0_xyy[k] * ab_z + g_zz_0_0_xyyz[k];

                g_zz_0_z_xyz[k] = -2.0 * g_z_0_0_xyz[k] - g_zz_0_0_xyz[k] * ab_z + g_zz_0_0_xyzz[k];

                g_zz_0_z_xzz[k] = -2.0 * g_z_0_0_xzz[k] - g_zz_0_0_xzz[k] * ab_z + g_zz_0_0_xzzz[k];

                g_zz_0_z_yyy[k] = -2.0 * g_z_0_0_yyy[k] - g_zz_0_0_yyy[k] * ab_z + g_zz_0_0_yyyz[k];

                g_zz_0_z_yyz[k] = -2.0 * g_z_0_0_yyz[k] - g_zz_0_0_yyz[k] * ab_z + g_zz_0_0_yyzz[k];

                g_zz_0_z_yzz[k] = -2.0 * g_z_0_0_yzz[k] - g_zz_0_0_yzz[k] * ab_z + g_zz_0_0_yzzz[k];

                g_zz_0_z_zzz[k] = -2.0 * g_z_0_0_zzz[k] - g_zz_0_0_zzz[k] * ab_z + g_zz_0_0_zzzz[k];
            }
        }
    }
}

} // erirec namespace

