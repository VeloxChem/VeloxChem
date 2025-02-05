#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_pfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_pfxx,
                                            const size_t idx_geom_01_sfxx,
                                            const size_t idx_geom_10_sfxx,
                                            const size_t idx_geom_11_sfxx,
                                            const size_t idx_geom_11_sgxx,
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

            const auto sf_geom_01_off = idx_geom_01_sfxx + i * dcomps + j;

            auto g_0_x_0_xxx = cbuffer.data(sf_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xxy = cbuffer.data(sf_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xxz = cbuffer.data(sf_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_xyy = cbuffer.data(sf_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_xyz = cbuffer.data(sf_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_xzz = cbuffer.data(sf_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_0_yyy = cbuffer.data(sf_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_0_yyz = cbuffer.data(sf_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_0_yzz = cbuffer.data(sf_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_0_zzz = cbuffer.data(sf_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_y_0_xxx = cbuffer.data(sf_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_y_0_xxy = cbuffer.data(sf_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_y_0_xxz = cbuffer.data(sf_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_y_0_xyy = cbuffer.data(sf_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_y_0_xyz = cbuffer.data(sf_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_y_0_xzz = cbuffer.data(sf_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_y_0_yyy = cbuffer.data(sf_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_y_0_yyz = cbuffer.data(sf_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_y_0_yzz = cbuffer.data(sf_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_0_zzz = cbuffer.data(sf_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_z_0_xxx = cbuffer.data(sf_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_z_0_xxy = cbuffer.data(sf_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_z_0_xxz = cbuffer.data(sf_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_z_0_xyy = cbuffer.data(sf_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_z_0_xyz = cbuffer.data(sf_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_z_0_xzz = cbuffer.data(sf_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_z_0_yyy = cbuffer.data(sf_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_z_0_yyz = cbuffer.data(sf_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_z_0_yzz = cbuffer.data(sf_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_z_0_zzz = cbuffer.data(sf_geom_01_off + 29 * ccomps * dcomps);

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

            const auto sf_geom_11_off = idx_geom_11_sfxx + i * dcomps + j;

            auto g_x_x_0_xxx = cbuffer.data(sf_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_0_xxy = cbuffer.data(sf_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_0_xxz = cbuffer.data(sf_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_0_xyy = cbuffer.data(sf_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_0_xyz = cbuffer.data(sf_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_0_xzz = cbuffer.data(sf_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_0_yyy = cbuffer.data(sf_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_0_yyz = cbuffer.data(sf_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_0_yzz = cbuffer.data(sf_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_0_zzz = cbuffer.data(sf_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_y_0_xxx = cbuffer.data(sf_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_y_0_xxy = cbuffer.data(sf_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_y_0_xxz = cbuffer.data(sf_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_y_0_xyy = cbuffer.data(sf_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_y_0_xyz = cbuffer.data(sf_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_y_0_xzz = cbuffer.data(sf_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_y_0_yyy = cbuffer.data(sf_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_y_0_yyz = cbuffer.data(sf_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_y_0_yzz = cbuffer.data(sf_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_y_0_zzz = cbuffer.data(sf_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_z_0_xxx = cbuffer.data(sf_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_z_0_xxy = cbuffer.data(sf_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_z_0_xxz = cbuffer.data(sf_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_z_0_xyy = cbuffer.data(sf_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_z_0_xyz = cbuffer.data(sf_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_z_0_xzz = cbuffer.data(sf_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_z_0_yyy = cbuffer.data(sf_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_z_0_yyz = cbuffer.data(sf_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_z_0_yzz = cbuffer.data(sf_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_z_0_zzz = cbuffer.data(sf_geom_11_off + 29 * ccomps * dcomps);

            auto g_y_x_0_xxx = cbuffer.data(sf_geom_11_off + 30 * ccomps * dcomps);

            auto g_y_x_0_xxy = cbuffer.data(sf_geom_11_off + 31 * ccomps * dcomps);

            auto g_y_x_0_xxz = cbuffer.data(sf_geom_11_off + 32 * ccomps * dcomps);

            auto g_y_x_0_xyy = cbuffer.data(sf_geom_11_off + 33 * ccomps * dcomps);

            auto g_y_x_0_xyz = cbuffer.data(sf_geom_11_off + 34 * ccomps * dcomps);

            auto g_y_x_0_xzz = cbuffer.data(sf_geom_11_off + 35 * ccomps * dcomps);

            auto g_y_x_0_yyy = cbuffer.data(sf_geom_11_off + 36 * ccomps * dcomps);

            auto g_y_x_0_yyz = cbuffer.data(sf_geom_11_off + 37 * ccomps * dcomps);

            auto g_y_x_0_yzz = cbuffer.data(sf_geom_11_off + 38 * ccomps * dcomps);

            auto g_y_x_0_zzz = cbuffer.data(sf_geom_11_off + 39 * ccomps * dcomps);

            auto g_y_y_0_xxx = cbuffer.data(sf_geom_11_off + 40 * ccomps * dcomps);

            auto g_y_y_0_xxy = cbuffer.data(sf_geom_11_off + 41 * ccomps * dcomps);

            auto g_y_y_0_xxz = cbuffer.data(sf_geom_11_off + 42 * ccomps * dcomps);

            auto g_y_y_0_xyy = cbuffer.data(sf_geom_11_off + 43 * ccomps * dcomps);

            auto g_y_y_0_xyz = cbuffer.data(sf_geom_11_off + 44 * ccomps * dcomps);

            auto g_y_y_0_xzz = cbuffer.data(sf_geom_11_off + 45 * ccomps * dcomps);

            auto g_y_y_0_yyy = cbuffer.data(sf_geom_11_off + 46 * ccomps * dcomps);

            auto g_y_y_0_yyz = cbuffer.data(sf_geom_11_off + 47 * ccomps * dcomps);

            auto g_y_y_0_yzz = cbuffer.data(sf_geom_11_off + 48 * ccomps * dcomps);

            auto g_y_y_0_zzz = cbuffer.data(sf_geom_11_off + 49 * ccomps * dcomps);

            auto g_y_z_0_xxx = cbuffer.data(sf_geom_11_off + 50 * ccomps * dcomps);

            auto g_y_z_0_xxy = cbuffer.data(sf_geom_11_off + 51 * ccomps * dcomps);

            auto g_y_z_0_xxz = cbuffer.data(sf_geom_11_off + 52 * ccomps * dcomps);

            auto g_y_z_0_xyy = cbuffer.data(sf_geom_11_off + 53 * ccomps * dcomps);

            auto g_y_z_0_xyz = cbuffer.data(sf_geom_11_off + 54 * ccomps * dcomps);

            auto g_y_z_0_xzz = cbuffer.data(sf_geom_11_off + 55 * ccomps * dcomps);

            auto g_y_z_0_yyy = cbuffer.data(sf_geom_11_off + 56 * ccomps * dcomps);

            auto g_y_z_0_yyz = cbuffer.data(sf_geom_11_off + 57 * ccomps * dcomps);

            auto g_y_z_0_yzz = cbuffer.data(sf_geom_11_off + 58 * ccomps * dcomps);

            auto g_y_z_0_zzz = cbuffer.data(sf_geom_11_off + 59 * ccomps * dcomps);

            auto g_z_x_0_xxx = cbuffer.data(sf_geom_11_off + 60 * ccomps * dcomps);

            auto g_z_x_0_xxy = cbuffer.data(sf_geom_11_off + 61 * ccomps * dcomps);

            auto g_z_x_0_xxz = cbuffer.data(sf_geom_11_off + 62 * ccomps * dcomps);

            auto g_z_x_0_xyy = cbuffer.data(sf_geom_11_off + 63 * ccomps * dcomps);

            auto g_z_x_0_xyz = cbuffer.data(sf_geom_11_off + 64 * ccomps * dcomps);

            auto g_z_x_0_xzz = cbuffer.data(sf_geom_11_off + 65 * ccomps * dcomps);

            auto g_z_x_0_yyy = cbuffer.data(sf_geom_11_off + 66 * ccomps * dcomps);

            auto g_z_x_0_yyz = cbuffer.data(sf_geom_11_off + 67 * ccomps * dcomps);

            auto g_z_x_0_yzz = cbuffer.data(sf_geom_11_off + 68 * ccomps * dcomps);

            auto g_z_x_0_zzz = cbuffer.data(sf_geom_11_off + 69 * ccomps * dcomps);

            auto g_z_y_0_xxx = cbuffer.data(sf_geom_11_off + 70 * ccomps * dcomps);

            auto g_z_y_0_xxy = cbuffer.data(sf_geom_11_off + 71 * ccomps * dcomps);

            auto g_z_y_0_xxz = cbuffer.data(sf_geom_11_off + 72 * ccomps * dcomps);

            auto g_z_y_0_xyy = cbuffer.data(sf_geom_11_off + 73 * ccomps * dcomps);

            auto g_z_y_0_xyz = cbuffer.data(sf_geom_11_off + 74 * ccomps * dcomps);

            auto g_z_y_0_xzz = cbuffer.data(sf_geom_11_off + 75 * ccomps * dcomps);

            auto g_z_y_0_yyy = cbuffer.data(sf_geom_11_off + 76 * ccomps * dcomps);

            auto g_z_y_0_yyz = cbuffer.data(sf_geom_11_off + 77 * ccomps * dcomps);

            auto g_z_y_0_yzz = cbuffer.data(sf_geom_11_off + 78 * ccomps * dcomps);

            auto g_z_y_0_zzz = cbuffer.data(sf_geom_11_off + 79 * ccomps * dcomps);

            auto g_z_z_0_xxx = cbuffer.data(sf_geom_11_off + 80 * ccomps * dcomps);

            auto g_z_z_0_xxy = cbuffer.data(sf_geom_11_off + 81 * ccomps * dcomps);

            auto g_z_z_0_xxz = cbuffer.data(sf_geom_11_off + 82 * ccomps * dcomps);

            auto g_z_z_0_xyy = cbuffer.data(sf_geom_11_off + 83 * ccomps * dcomps);

            auto g_z_z_0_xyz = cbuffer.data(sf_geom_11_off + 84 * ccomps * dcomps);

            auto g_z_z_0_xzz = cbuffer.data(sf_geom_11_off + 85 * ccomps * dcomps);

            auto g_z_z_0_yyy = cbuffer.data(sf_geom_11_off + 86 * ccomps * dcomps);

            auto g_z_z_0_yyz = cbuffer.data(sf_geom_11_off + 87 * ccomps * dcomps);

            auto g_z_z_0_yzz = cbuffer.data(sf_geom_11_off + 88 * ccomps * dcomps);

            auto g_z_z_0_zzz = cbuffer.data(sf_geom_11_off + 89 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SGSS

            const auto sg_geom_11_off = idx_geom_11_sgxx + i * dcomps + j;

            auto g_x_x_0_xxxx = cbuffer.data(sg_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_0_xxxy = cbuffer.data(sg_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_0_xxxz = cbuffer.data(sg_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_0_xxyy = cbuffer.data(sg_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_0_xxyz = cbuffer.data(sg_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_0_xxzz = cbuffer.data(sg_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_0_xyyy = cbuffer.data(sg_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_0_xyyz = cbuffer.data(sg_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_0_xyzz = cbuffer.data(sg_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_0_xzzz = cbuffer.data(sg_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_0_yyyy = cbuffer.data(sg_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_0_yyyz = cbuffer.data(sg_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_0_yyzz = cbuffer.data(sg_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_0_yzzz = cbuffer.data(sg_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_0_zzzz = cbuffer.data(sg_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_y_0_xxxx = cbuffer.data(sg_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_y_0_xxxy = cbuffer.data(sg_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_y_0_xxxz = cbuffer.data(sg_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_y_0_xxyy = cbuffer.data(sg_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_y_0_xxyz = cbuffer.data(sg_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_y_0_xxzz = cbuffer.data(sg_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_y_0_xyyy = cbuffer.data(sg_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_y_0_xyyz = cbuffer.data(sg_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_y_0_xyzz = cbuffer.data(sg_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_y_0_xzzz = cbuffer.data(sg_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_y_0_yyyy = cbuffer.data(sg_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_y_0_yyyz = cbuffer.data(sg_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_y_0_yyzz = cbuffer.data(sg_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_y_0_yzzz = cbuffer.data(sg_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_y_0_zzzz = cbuffer.data(sg_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_z_0_xxxx = cbuffer.data(sg_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_z_0_xxxy = cbuffer.data(sg_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_z_0_xxxz = cbuffer.data(sg_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_z_0_xxyy = cbuffer.data(sg_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_z_0_xxyz = cbuffer.data(sg_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_z_0_xxzz = cbuffer.data(sg_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_z_0_xyyy = cbuffer.data(sg_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_z_0_xyyz = cbuffer.data(sg_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_z_0_xyzz = cbuffer.data(sg_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_z_0_xzzz = cbuffer.data(sg_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_z_0_yyyy = cbuffer.data(sg_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_z_0_yyyz = cbuffer.data(sg_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_z_0_yyzz = cbuffer.data(sg_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_z_0_yzzz = cbuffer.data(sg_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_z_0_zzzz = cbuffer.data(sg_geom_11_off + 44 * ccomps * dcomps);

            auto g_y_x_0_xxxx = cbuffer.data(sg_geom_11_off + 45 * ccomps * dcomps);

            auto g_y_x_0_xxxy = cbuffer.data(sg_geom_11_off + 46 * ccomps * dcomps);

            auto g_y_x_0_xxxz = cbuffer.data(sg_geom_11_off + 47 * ccomps * dcomps);

            auto g_y_x_0_xxyy = cbuffer.data(sg_geom_11_off + 48 * ccomps * dcomps);

            auto g_y_x_0_xxyz = cbuffer.data(sg_geom_11_off + 49 * ccomps * dcomps);

            auto g_y_x_0_xxzz = cbuffer.data(sg_geom_11_off + 50 * ccomps * dcomps);

            auto g_y_x_0_xyyy = cbuffer.data(sg_geom_11_off + 51 * ccomps * dcomps);

            auto g_y_x_0_xyyz = cbuffer.data(sg_geom_11_off + 52 * ccomps * dcomps);

            auto g_y_x_0_xyzz = cbuffer.data(sg_geom_11_off + 53 * ccomps * dcomps);

            auto g_y_x_0_xzzz = cbuffer.data(sg_geom_11_off + 54 * ccomps * dcomps);

            auto g_y_x_0_yyyy = cbuffer.data(sg_geom_11_off + 55 * ccomps * dcomps);

            auto g_y_x_0_yyyz = cbuffer.data(sg_geom_11_off + 56 * ccomps * dcomps);

            auto g_y_x_0_yyzz = cbuffer.data(sg_geom_11_off + 57 * ccomps * dcomps);

            auto g_y_x_0_yzzz = cbuffer.data(sg_geom_11_off + 58 * ccomps * dcomps);

            auto g_y_x_0_zzzz = cbuffer.data(sg_geom_11_off + 59 * ccomps * dcomps);

            auto g_y_y_0_xxxx = cbuffer.data(sg_geom_11_off + 60 * ccomps * dcomps);

            auto g_y_y_0_xxxy = cbuffer.data(sg_geom_11_off + 61 * ccomps * dcomps);

            auto g_y_y_0_xxxz = cbuffer.data(sg_geom_11_off + 62 * ccomps * dcomps);

            auto g_y_y_0_xxyy = cbuffer.data(sg_geom_11_off + 63 * ccomps * dcomps);

            auto g_y_y_0_xxyz = cbuffer.data(sg_geom_11_off + 64 * ccomps * dcomps);

            auto g_y_y_0_xxzz = cbuffer.data(sg_geom_11_off + 65 * ccomps * dcomps);

            auto g_y_y_0_xyyy = cbuffer.data(sg_geom_11_off + 66 * ccomps * dcomps);

            auto g_y_y_0_xyyz = cbuffer.data(sg_geom_11_off + 67 * ccomps * dcomps);

            auto g_y_y_0_xyzz = cbuffer.data(sg_geom_11_off + 68 * ccomps * dcomps);

            auto g_y_y_0_xzzz = cbuffer.data(sg_geom_11_off + 69 * ccomps * dcomps);

            auto g_y_y_0_yyyy = cbuffer.data(sg_geom_11_off + 70 * ccomps * dcomps);

            auto g_y_y_0_yyyz = cbuffer.data(sg_geom_11_off + 71 * ccomps * dcomps);

            auto g_y_y_0_yyzz = cbuffer.data(sg_geom_11_off + 72 * ccomps * dcomps);

            auto g_y_y_0_yzzz = cbuffer.data(sg_geom_11_off + 73 * ccomps * dcomps);

            auto g_y_y_0_zzzz = cbuffer.data(sg_geom_11_off + 74 * ccomps * dcomps);

            auto g_y_z_0_xxxx = cbuffer.data(sg_geom_11_off + 75 * ccomps * dcomps);

            auto g_y_z_0_xxxy = cbuffer.data(sg_geom_11_off + 76 * ccomps * dcomps);

            auto g_y_z_0_xxxz = cbuffer.data(sg_geom_11_off + 77 * ccomps * dcomps);

            auto g_y_z_0_xxyy = cbuffer.data(sg_geom_11_off + 78 * ccomps * dcomps);

            auto g_y_z_0_xxyz = cbuffer.data(sg_geom_11_off + 79 * ccomps * dcomps);

            auto g_y_z_0_xxzz = cbuffer.data(sg_geom_11_off + 80 * ccomps * dcomps);

            auto g_y_z_0_xyyy = cbuffer.data(sg_geom_11_off + 81 * ccomps * dcomps);

            auto g_y_z_0_xyyz = cbuffer.data(sg_geom_11_off + 82 * ccomps * dcomps);

            auto g_y_z_0_xyzz = cbuffer.data(sg_geom_11_off + 83 * ccomps * dcomps);

            auto g_y_z_0_xzzz = cbuffer.data(sg_geom_11_off + 84 * ccomps * dcomps);

            auto g_y_z_0_yyyy = cbuffer.data(sg_geom_11_off + 85 * ccomps * dcomps);

            auto g_y_z_0_yyyz = cbuffer.data(sg_geom_11_off + 86 * ccomps * dcomps);

            auto g_y_z_0_yyzz = cbuffer.data(sg_geom_11_off + 87 * ccomps * dcomps);

            auto g_y_z_0_yzzz = cbuffer.data(sg_geom_11_off + 88 * ccomps * dcomps);

            auto g_y_z_0_zzzz = cbuffer.data(sg_geom_11_off + 89 * ccomps * dcomps);

            auto g_z_x_0_xxxx = cbuffer.data(sg_geom_11_off + 90 * ccomps * dcomps);

            auto g_z_x_0_xxxy = cbuffer.data(sg_geom_11_off + 91 * ccomps * dcomps);

            auto g_z_x_0_xxxz = cbuffer.data(sg_geom_11_off + 92 * ccomps * dcomps);

            auto g_z_x_0_xxyy = cbuffer.data(sg_geom_11_off + 93 * ccomps * dcomps);

            auto g_z_x_0_xxyz = cbuffer.data(sg_geom_11_off + 94 * ccomps * dcomps);

            auto g_z_x_0_xxzz = cbuffer.data(sg_geom_11_off + 95 * ccomps * dcomps);

            auto g_z_x_0_xyyy = cbuffer.data(sg_geom_11_off + 96 * ccomps * dcomps);

            auto g_z_x_0_xyyz = cbuffer.data(sg_geom_11_off + 97 * ccomps * dcomps);

            auto g_z_x_0_xyzz = cbuffer.data(sg_geom_11_off + 98 * ccomps * dcomps);

            auto g_z_x_0_xzzz = cbuffer.data(sg_geom_11_off + 99 * ccomps * dcomps);

            auto g_z_x_0_yyyy = cbuffer.data(sg_geom_11_off + 100 * ccomps * dcomps);

            auto g_z_x_0_yyyz = cbuffer.data(sg_geom_11_off + 101 * ccomps * dcomps);

            auto g_z_x_0_yyzz = cbuffer.data(sg_geom_11_off + 102 * ccomps * dcomps);

            auto g_z_x_0_yzzz = cbuffer.data(sg_geom_11_off + 103 * ccomps * dcomps);

            auto g_z_x_0_zzzz = cbuffer.data(sg_geom_11_off + 104 * ccomps * dcomps);

            auto g_z_y_0_xxxx = cbuffer.data(sg_geom_11_off + 105 * ccomps * dcomps);

            auto g_z_y_0_xxxy = cbuffer.data(sg_geom_11_off + 106 * ccomps * dcomps);

            auto g_z_y_0_xxxz = cbuffer.data(sg_geom_11_off + 107 * ccomps * dcomps);

            auto g_z_y_0_xxyy = cbuffer.data(sg_geom_11_off + 108 * ccomps * dcomps);

            auto g_z_y_0_xxyz = cbuffer.data(sg_geom_11_off + 109 * ccomps * dcomps);

            auto g_z_y_0_xxzz = cbuffer.data(sg_geom_11_off + 110 * ccomps * dcomps);

            auto g_z_y_0_xyyy = cbuffer.data(sg_geom_11_off + 111 * ccomps * dcomps);

            auto g_z_y_0_xyyz = cbuffer.data(sg_geom_11_off + 112 * ccomps * dcomps);

            auto g_z_y_0_xyzz = cbuffer.data(sg_geom_11_off + 113 * ccomps * dcomps);

            auto g_z_y_0_xzzz = cbuffer.data(sg_geom_11_off + 114 * ccomps * dcomps);

            auto g_z_y_0_yyyy = cbuffer.data(sg_geom_11_off + 115 * ccomps * dcomps);

            auto g_z_y_0_yyyz = cbuffer.data(sg_geom_11_off + 116 * ccomps * dcomps);

            auto g_z_y_0_yyzz = cbuffer.data(sg_geom_11_off + 117 * ccomps * dcomps);

            auto g_z_y_0_yzzz = cbuffer.data(sg_geom_11_off + 118 * ccomps * dcomps);

            auto g_z_y_0_zzzz = cbuffer.data(sg_geom_11_off + 119 * ccomps * dcomps);

            auto g_z_z_0_xxxx = cbuffer.data(sg_geom_11_off + 120 * ccomps * dcomps);

            auto g_z_z_0_xxxy = cbuffer.data(sg_geom_11_off + 121 * ccomps * dcomps);

            auto g_z_z_0_xxxz = cbuffer.data(sg_geom_11_off + 122 * ccomps * dcomps);

            auto g_z_z_0_xxyy = cbuffer.data(sg_geom_11_off + 123 * ccomps * dcomps);

            auto g_z_z_0_xxyz = cbuffer.data(sg_geom_11_off + 124 * ccomps * dcomps);

            auto g_z_z_0_xxzz = cbuffer.data(sg_geom_11_off + 125 * ccomps * dcomps);

            auto g_z_z_0_xyyy = cbuffer.data(sg_geom_11_off + 126 * ccomps * dcomps);

            auto g_z_z_0_xyyz = cbuffer.data(sg_geom_11_off + 127 * ccomps * dcomps);

            auto g_z_z_0_xyzz = cbuffer.data(sg_geom_11_off + 128 * ccomps * dcomps);

            auto g_z_z_0_xzzz = cbuffer.data(sg_geom_11_off + 129 * ccomps * dcomps);

            auto g_z_z_0_yyyy = cbuffer.data(sg_geom_11_off + 130 * ccomps * dcomps);

            auto g_z_z_0_yyyz = cbuffer.data(sg_geom_11_off + 131 * ccomps * dcomps);

            auto g_z_z_0_yyzz = cbuffer.data(sg_geom_11_off + 132 * ccomps * dcomps);

            auto g_z_z_0_yzzz = cbuffer.data(sg_geom_11_off + 133 * ccomps * dcomps);

            auto g_z_z_0_zzzz = cbuffer.data(sg_geom_11_off + 134 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pfxx

            const auto pf_geom_11_off = idx_geom_11_pfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_x_x_xxx = cbuffer.data(pf_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xxy = cbuffer.data(pf_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xxz = cbuffer.data(pf_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_xyy = cbuffer.data(pf_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_xyz = cbuffer.data(pf_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_xzz = cbuffer.data(pf_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_x_yyy = cbuffer.data(pf_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_x_yyz = cbuffer.data(pf_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_x_yzz = cbuffer.data(pf_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_x_zzz = cbuffer.data(pf_geom_11_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxx, g_0_x_0_xxy, g_0_x_0_xxz, g_0_x_0_xyy, g_0_x_0_xyz, g_0_x_0_xzz, g_0_x_0_yyy, g_0_x_0_yyz, g_0_x_0_yzz, g_0_x_0_zzz, g_x_0_0_xxx, g_x_0_0_xxy, g_x_0_0_xxz, g_x_0_0_xyy, g_x_0_0_xyz, g_x_0_0_xzz, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yzz, g_x_0_0_zzz, g_x_x_0_xxx, g_x_x_0_xxxx, g_x_x_0_xxxy, g_x_x_0_xxxz, g_x_x_0_xxy, g_x_x_0_xxyy, g_x_x_0_xxyz, g_x_x_0_xxz, g_x_x_0_xxzz, g_x_x_0_xyy, g_x_x_0_xyyy, g_x_x_0_xyyz, g_x_x_0_xyz, g_x_x_0_xyzz, g_x_x_0_xzz, g_x_x_0_xzzz, g_x_x_0_yyy, g_x_x_0_yyz, g_x_x_0_yzz, g_x_x_0_zzz, g_x_x_x_xxx, g_x_x_x_xxy, g_x_x_x_xxz, g_x_x_x_xyy, g_x_x_x_xyz, g_x_x_x_xzz, g_x_x_x_yyy, g_x_x_x_yyz, g_x_x_x_yzz, g_x_x_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_x_xxx[k] = -g_0_x_0_xxx[k] + g_x_0_0_xxx[k] - g_x_x_0_xxx[k] * ab_x + g_x_x_0_xxxx[k];

                g_x_x_x_xxy[k] = -g_0_x_0_xxy[k] + g_x_0_0_xxy[k] - g_x_x_0_xxy[k] * ab_x + g_x_x_0_xxxy[k];

                g_x_x_x_xxz[k] = -g_0_x_0_xxz[k] + g_x_0_0_xxz[k] - g_x_x_0_xxz[k] * ab_x + g_x_x_0_xxxz[k];

                g_x_x_x_xyy[k] = -g_0_x_0_xyy[k] + g_x_0_0_xyy[k] - g_x_x_0_xyy[k] * ab_x + g_x_x_0_xxyy[k];

                g_x_x_x_xyz[k] = -g_0_x_0_xyz[k] + g_x_0_0_xyz[k] - g_x_x_0_xyz[k] * ab_x + g_x_x_0_xxyz[k];

                g_x_x_x_xzz[k] = -g_0_x_0_xzz[k] + g_x_0_0_xzz[k] - g_x_x_0_xzz[k] * ab_x + g_x_x_0_xxzz[k];

                g_x_x_x_yyy[k] = -g_0_x_0_yyy[k] + g_x_0_0_yyy[k] - g_x_x_0_yyy[k] * ab_x + g_x_x_0_xyyy[k];

                g_x_x_x_yyz[k] = -g_0_x_0_yyz[k] + g_x_0_0_yyz[k] - g_x_x_0_yyz[k] * ab_x + g_x_x_0_xyyz[k];

                g_x_x_x_yzz[k] = -g_0_x_0_yzz[k] + g_x_0_0_yzz[k] - g_x_x_0_yzz[k] * ab_x + g_x_x_0_xyzz[k];

                g_x_x_x_zzz[k] = -g_0_x_0_zzz[k] + g_x_0_0_zzz[k] - g_x_x_0_zzz[k] * ab_x + g_x_x_0_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_x_y_xxx = cbuffer.data(pf_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_y_xxy = cbuffer.data(pf_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_y_xxz = cbuffer.data(pf_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_y_xyy = cbuffer.data(pf_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_y_xyz = cbuffer.data(pf_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_y_xzz = cbuffer.data(pf_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_y_yyy = cbuffer.data(pf_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_y_yyz = cbuffer.data(pf_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_y_yzz = cbuffer.data(pf_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_y_zzz = cbuffer.data(pf_geom_11_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_0_xxx, g_x_x_0_xxxy, g_x_x_0_xxy, g_x_x_0_xxyy, g_x_x_0_xxyz, g_x_x_0_xxz, g_x_x_0_xyy, g_x_x_0_xyyy, g_x_x_0_xyyz, g_x_x_0_xyz, g_x_x_0_xyzz, g_x_x_0_xzz, g_x_x_0_yyy, g_x_x_0_yyyy, g_x_x_0_yyyz, g_x_x_0_yyz, g_x_x_0_yyzz, g_x_x_0_yzz, g_x_x_0_yzzz, g_x_x_0_zzz, g_x_x_y_xxx, g_x_x_y_xxy, g_x_x_y_xxz, g_x_x_y_xyy, g_x_x_y_xyz, g_x_x_y_xzz, g_x_x_y_yyy, g_x_x_y_yyz, g_x_x_y_yzz, g_x_x_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_y_xxx[k] = -g_x_x_0_xxx[k] * ab_y + g_x_x_0_xxxy[k];

                g_x_x_y_xxy[k] = -g_x_x_0_xxy[k] * ab_y + g_x_x_0_xxyy[k];

                g_x_x_y_xxz[k] = -g_x_x_0_xxz[k] * ab_y + g_x_x_0_xxyz[k];

                g_x_x_y_xyy[k] = -g_x_x_0_xyy[k] * ab_y + g_x_x_0_xyyy[k];

                g_x_x_y_xyz[k] = -g_x_x_0_xyz[k] * ab_y + g_x_x_0_xyyz[k];

                g_x_x_y_xzz[k] = -g_x_x_0_xzz[k] * ab_y + g_x_x_0_xyzz[k];

                g_x_x_y_yyy[k] = -g_x_x_0_yyy[k] * ab_y + g_x_x_0_yyyy[k];

                g_x_x_y_yyz[k] = -g_x_x_0_yyz[k] * ab_y + g_x_x_0_yyyz[k];

                g_x_x_y_yzz[k] = -g_x_x_0_yzz[k] * ab_y + g_x_x_0_yyzz[k];

                g_x_x_y_zzz[k] = -g_x_x_0_zzz[k] * ab_y + g_x_x_0_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_x_z_xxx = cbuffer.data(pf_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_z_xxy = cbuffer.data(pf_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_z_xxz = cbuffer.data(pf_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_z_xyy = cbuffer.data(pf_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_z_xyz = cbuffer.data(pf_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_z_xzz = cbuffer.data(pf_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_z_yyy = cbuffer.data(pf_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_z_yyz = cbuffer.data(pf_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_z_yzz = cbuffer.data(pf_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_z_zzz = cbuffer.data(pf_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_0_xxx, g_x_x_0_xxxz, g_x_x_0_xxy, g_x_x_0_xxyz, g_x_x_0_xxz, g_x_x_0_xxzz, g_x_x_0_xyy, g_x_x_0_xyyz, g_x_x_0_xyz, g_x_x_0_xyzz, g_x_x_0_xzz, g_x_x_0_xzzz, g_x_x_0_yyy, g_x_x_0_yyyz, g_x_x_0_yyz, g_x_x_0_yyzz, g_x_x_0_yzz, g_x_x_0_yzzz, g_x_x_0_zzz, g_x_x_0_zzzz, g_x_x_z_xxx, g_x_x_z_xxy, g_x_x_z_xxz, g_x_x_z_xyy, g_x_x_z_xyz, g_x_x_z_xzz, g_x_x_z_yyy, g_x_x_z_yyz, g_x_x_z_yzz, g_x_x_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_z_xxx[k] = -g_x_x_0_xxx[k] * ab_z + g_x_x_0_xxxz[k];

                g_x_x_z_xxy[k] = -g_x_x_0_xxy[k] * ab_z + g_x_x_0_xxyz[k];

                g_x_x_z_xxz[k] = -g_x_x_0_xxz[k] * ab_z + g_x_x_0_xxzz[k];

                g_x_x_z_xyy[k] = -g_x_x_0_xyy[k] * ab_z + g_x_x_0_xyyz[k];

                g_x_x_z_xyz[k] = -g_x_x_0_xyz[k] * ab_z + g_x_x_0_xyzz[k];

                g_x_x_z_xzz[k] = -g_x_x_0_xzz[k] * ab_z + g_x_x_0_xzzz[k];

                g_x_x_z_yyy[k] = -g_x_x_0_yyy[k] * ab_z + g_x_x_0_yyyz[k];

                g_x_x_z_yyz[k] = -g_x_x_0_yyz[k] * ab_z + g_x_x_0_yyzz[k];

                g_x_x_z_yzz[k] = -g_x_x_0_yzz[k] * ab_z + g_x_x_0_yzzz[k];

                g_x_x_z_zzz[k] = -g_x_x_0_zzz[k] * ab_z + g_x_x_0_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_x_y_x_xxx = cbuffer.data(pf_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_y_x_xxy = cbuffer.data(pf_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_y_x_xxz = cbuffer.data(pf_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_y_x_xyy = cbuffer.data(pf_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_y_x_xyz = cbuffer.data(pf_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_x_xzz = cbuffer.data(pf_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_y_x_yyy = cbuffer.data(pf_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_y_x_yyz = cbuffer.data(pf_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_y_x_yzz = cbuffer.data(pf_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_y_x_zzz = cbuffer.data(pf_geom_11_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxx, g_0_y_0_xxy, g_0_y_0_xxz, g_0_y_0_xyy, g_0_y_0_xyz, g_0_y_0_xzz, g_0_y_0_yyy, g_0_y_0_yyz, g_0_y_0_yzz, g_0_y_0_zzz, g_x_y_0_xxx, g_x_y_0_xxxx, g_x_y_0_xxxy, g_x_y_0_xxxz, g_x_y_0_xxy, g_x_y_0_xxyy, g_x_y_0_xxyz, g_x_y_0_xxz, g_x_y_0_xxzz, g_x_y_0_xyy, g_x_y_0_xyyy, g_x_y_0_xyyz, g_x_y_0_xyz, g_x_y_0_xyzz, g_x_y_0_xzz, g_x_y_0_xzzz, g_x_y_0_yyy, g_x_y_0_yyz, g_x_y_0_yzz, g_x_y_0_zzz, g_x_y_x_xxx, g_x_y_x_xxy, g_x_y_x_xxz, g_x_y_x_xyy, g_x_y_x_xyz, g_x_y_x_xzz, g_x_y_x_yyy, g_x_y_x_yyz, g_x_y_x_yzz, g_x_y_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_x_xxx[k] = -g_0_y_0_xxx[k] - g_x_y_0_xxx[k] * ab_x + g_x_y_0_xxxx[k];

                g_x_y_x_xxy[k] = -g_0_y_0_xxy[k] - g_x_y_0_xxy[k] * ab_x + g_x_y_0_xxxy[k];

                g_x_y_x_xxz[k] = -g_0_y_0_xxz[k] - g_x_y_0_xxz[k] * ab_x + g_x_y_0_xxxz[k];

                g_x_y_x_xyy[k] = -g_0_y_0_xyy[k] - g_x_y_0_xyy[k] * ab_x + g_x_y_0_xxyy[k];

                g_x_y_x_xyz[k] = -g_0_y_0_xyz[k] - g_x_y_0_xyz[k] * ab_x + g_x_y_0_xxyz[k];

                g_x_y_x_xzz[k] = -g_0_y_0_xzz[k] - g_x_y_0_xzz[k] * ab_x + g_x_y_0_xxzz[k];

                g_x_y_x_yyy[k] = -g_0_y_0_yyy[k] - g_x_y_0_yyy[k] * ab_x + g_x_y_0_xyyy[k];

                g_x_y_x_yyz[k] = -g_0_y_0_yyz[k] - g_x_y_0_yyz[k] * ab_x + g_x_y_0_xyyz[k];

                g_x_y_x_yzz[k] = -g_0_y_0_yzz[k] - g_x_y_0_yzz[k] * ab_x + g_x_y_0_xyzz[k];

                g_x_y_x_zzz[k] = -g_0_y_0_zzz[k] - g_x_y_0_zzz[k] * ab_x + g_x_y_0_xzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_x_y_y_xxx = cbuffer.data(pf_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_y_y_xxy = cbuffer.data(pf_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_y_y_xxz = cbuffer.data(pf_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_y_y_xyy = cbuffer.data(pf_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_y_y_xyz = cbuffer.data(pf_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_y_xzz = cbuffer.data(pf_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_y_yyy = cbuffer.data(pf_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_y_yyz = cbuffer.data(pf_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_y_y_yzz = cbuffer.data(pf_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_y_zzz = cbuffer.data(pf_geom_11_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxx, g_x_0_0_xxy, g_x_0_0_xxz, g_x_0_0_xyy, g_x_0_0_xyz, g_x_0_0_xzz, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yzz, g_x_0_0_zzz, g_x_y_0_xxx, g_x_y_0_xxxy, g_x_y_0_xxy, g_x_y_0_xxyy, g_x_y_0_xxyz, g_x_y_0_xxz, g_x_y_0_xyy, g_x_y_0_xyyy, g_x_y_0_xyyz, g_x_y_0_xyz, g_x_y_0_xyzz, g_x_y_0_xzz, g_x_y_0_yyy, g_x_y_0_yyyy, g_x_y_0_yyyz, g_x_y_0_yyz, g_x_y_0_yyzz, g_x_y_0_yzz, g_x_y_0_yzzz, g_x_y_0_zzz, g_x_y_y_xxx, g_x_y_y_xxy, g_x_y_y_xxz, g_x_y_y_xyy, g_x_y_y_xyz, g_x_y_y_xzz, g_x_y_y_yyy, g_x_y_y_yyz, g_x_y_y_yzz, g_x_y_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_y_xxx[k] = g_x_0_0_xxx[k] - g_x_y_0_xxx[k] * ab_y + g_x_y_0_xxxy[k];

                g_x_y_y_xxy[k] = g_x_0_0_xxy[k] - g_x_y_0_xxy[k] * ab_y + g_x_y_0_xxyy[k];

                g_x_y_y_xxz[k] = g_x_0_0_xxz[k] - g_x_y_0_xxz[k] * ab_y + g_x_y_0_xxyz[k];

                g_x_y_y_xyy[k] = g_x_0_0_xyy[k] - g_x_y_0_xyy[k] * ab_y + g_x_y_0_xyyy[k];

                g_x_y_y_xyz[k] = g_x_0_0_xyz[k] - g_x_y_0_xyz[k] * ab_y + g_x_y_0_xyyz[k];

                g_x_y_y_xzz[k] = g_x_0_0_xzz[k] - g_x_y_0_xzz[k] * ab_y + g_x_y_0_xyzz[k];

                g_x_y_y_yyy[k] = g_x_0_0_yyy[k] - g_x_y_0_yyy[k] * ab_y + g_x_y_0_yyyy[k];

                g_x_y_y_yyz[k] = g_x_0_0_yyz[k] - g_x_y_0_yyz[k] * ab_y + g_x_y_0_yyyz[k];

                g_x_y_y_yzz[k] = g_x_0_0_yzz[k] - g_x_y_0_yzz[k] * ab_y + g_x_y_0_yyzz[k];

                g_x_y_y_zzz[k] = g_x_0_0_zzz[k] - g_x_y_0_zzz[k] * ab_y + g_x_y_0_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_x_y_z_xxx = cbuffer.data(pf_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_z_xxy = cbuffer.data(pf_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_z_xxz = cbuffer.data(pf_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_z_xyy = cbuffer.data(pf_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_y_z_xyz = cbuffer.data(pf_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_z_xzz = cbuffer.data(pf_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_z_yyy = cbuffer.data(pf_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_y_z_yyz = cbuffer.data(pf_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_z_yzz = cbuffer.data(pf_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_z_zzz = cbuffer.data(pf_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_0_xxx, g_x_y_0_xxxz, g_x_y_0_xxy, g_x_y_0_xxyz, g_x_y_0_xxz, g_x_y_0_xxzz, g_x_y_0_xyy, g_x_y_0_xyyz, g_x_y_0_xyz, g_x_y_0_xyzz, g_x_y_0_xzz, g_x_y_0_xzzz, g_x_y_0_yyy, g_x_y_0_yyyz, g_x_y_0_yyz, g_x_y_0_yyzz, g_x_y_0_yzz, g_x_y_0_yzzz, g_x_y_0_zzz, g_x_y_0_zzzz, g_x_y_z_xxx, g_x_y_z_xxy, g_x_y_z_xxz, g_x_y_z_xyy, g_x_y_z_xyz, g_x_y_z_xzz, g_x_y_z_yyy, g_x_y_z_yyz, g_x_y_z_yzz, g_x_y_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_z_xxx[k] = -g_x_y_0_xxx[k] * ab_z + g_x_y_0_xxxz[k];

                g_x_y_z_xxy[k] = -g_x_y_0_xxy[k] * ab_z + g_x_y_0_xxyz[k];

                g_x_y_z_xxz[k] = -g_x_y_0_xxz[k] * ab_z + g_x_y_0_xxzz[k];

                g_x_y_z_xyy[k] = -g_x_y_0_xyy[k] * ab_z + g_x_y_0_xyyz[k];

                g_x_y_z_xyz[k] = -g_x_y_0_xyz[k] * ab_z + g_x_y_0_xyzz[k];

                g_x_y_z_xzz[k] = -g_x_y_0_xzz[k] * ab_z + g_x_y_0_xzzz[k];

                g_x_y_z_yyy[k] = -g_x_y_0_yyy[k] * ab_z + g_x_y_0_yyyz[k];

                g_x_y_z_yyz[k] = -g_x_y_0_yyz[k] * ab_z + g_x_y_0_yyzz[k];

                g_x_y_z_yzz[k] = -g_x_y_0_yzz[k] * ab_z + g_x_y_0_yzzz[k];

                g_x_y_z_zzz[k] = -g_x_y_0_zzz[k] * ab_z + g_x_y_0_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_x_z_x_xxx = cbuffer.data(pf_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_z_x_xxy = cbuffer.data(pf_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_z_x_xxz = cbuffer.data(pf_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_z_x_xyy = cbuffer.data(pf_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_z_x_xyz = cbuffer.data(pf_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_z_x_xzz = cbuffer.data(pf_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_z_x_yyy = cbuffer.data(pf_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_z_x_yyz = cbuffer.data(pf_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_z_x_yzz = cbuffer.data(pf_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_z_x_zzz = cbuffer.data(pf_geom_11_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxx, g_0_z_0_xxy, g_0_z_0_xxz, g_0_z_0_xyy, g_0_z_0_xyz, g_0_z_0_xzz, g_0_z_0_yyy, g_0_z_0_yyz, g_0_z_0_yzz, g_0_z_0_zzz, g_x_z_0_xxx, g_x_z_0_xxxx, g_x_z_0_xxxy, g_x_z_0_xxxz, g_x_z_0_xxy, g_x_z_0_xxyy, g_x_z_0_xxyz, g_x_z_0_xxz, g_x_z_0_xxzz, g_x_z_0_xyy, g_x_z_0_xyyy, g_x_z_0_xyyz, g_x_z_0_xyz, g_x_z_0_xyzz, g_x_z_0_xzz, g_x_z_0_xzzz, g_x_z_0_yyy, g_x_z_0_yyz, g_x_z_0_yzz, g_x_z_0_zzz, g_x_z_x_xxx, g_x_z_x_xxy, g_x_z_x_xxz, g_x_z_x_xyy, g_x_z_x_xyz, g_x_z_x_xzz, g_x_z_x_yyy, g_x_z_x_yyz, g_x_z_x_yzz, g_x_z_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_x_xxx[k] = -g_0_z_0_xxx[k] - g_x_z_0_xxx[k] * ab_x + g_x_z_0_xxxx[k];

                g_x_z_x_xxy[k] = -g_0_z_0_xxy[k] - g_x_z_0_xxy[k] * ab_x + g_x_z_0_xxxy[k];

                g_x_z_x_xxz[k] = -g_0_z_0_xxz[k] - g_x_z_0_xxz[k] * ab_x + g_x_z_0_xxxz[k];

                g_x_z_x_xyy[k] = -g_0_z_0_xyy[k] - g_x_z_0_xyy[k] * ab_x + g_x_z_0_xxyy[k];

                g_x_z_x_xyz[k] = -g_0_z_0_xyz[k] - g_x_z_0_xyz[k] * ab_x + g_x_z_0_xxyz[k];

                g_x_z_x_xzz[k] = -g_0_z_0_xzz[k] - g_x_z_0_xzz[k] * ab_x + g_x_z_0_xxzz[k];

                g_x_z_x_yyy[k] = -g_0_z_0_yyy[k] - g_x_z_0_yyy[k] * ab_x + g_x_z_0_xyyy[k];

                g_x_z_x_yyz[k] = -g_0_z_0_yyz[k] - g_x_z_0_yyz[k] * ab_x + g_x_z_0_xyyz[k];

                g_x_z_x_yzz[k] = -g_0_z_0_yzz[k] - g_x_z_0_yzz[k] * ab_x + g_x_z_0_xyzz[k];

                g_x_z_x_zzz[k] = -g_0_z_0_zzz[k] - g_x_z_0_zzz[k] * ab_x + g_x_z_0_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_x_z_y_xxx = cbuffer.data(pf_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_z_y_xxy = cbuffer.data(pf_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_z_y_xxz = cbuffer.data(pf_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_z_y_xyy = cbuffer.data(pf_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_z_y_xyz = cbuffer.data(pf_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_z_y_xzz = cbuffer.data(pf_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_z_y_yyy = cbuffer.data(pf_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_z_y_yyz = cbuffer.data(pf_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_z_y_yzz = cbuffer.data(pf_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_z_y_zzz = cbuffer.data(pf_geom_11_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_0_xxx, g_x_z_0_xxxy, g_x_z_0_xxy, g_x_z_0_xxyy, g_x_z_0_xxyz, g_x_z_0_xxz, g_x_z_0_xyy, g_x_z_0_xyyy, g_x_z_0_xyyz, g_x_z_0_xyz, g_x_z_0_xyzz, g_x_z_0_xzz, g_x_z_0_yyy, g_x_z_0_yyyy, g_x_z_0_yyyz, g_x_z_0_yyz, g_x_z_0_yyzz, g_x_z_0_yzz, g_x_z_0_yzzz, g_x_z_0_zzz, g_x_z_y_xxx, g_x_z_y_xxy, g_x_z_y_xxz, g_x_z_y_xyy, g_x_z_y_xyz, g_x_z_y_xzz, g_x_z_y_yyy, g_x_z_y_yyz, g_x_z_y_yzz, g_x_z_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_y_xxx[k] = -g_x_z_0_xxx[k] * ab_y + g_x_z_0_xxxy[k];

                g_x_z_y_xxy[k] = -g_x_z_0_xxy[k] * ab_y + g_x_z_0_xxyy[k];

                g_x_z_y_xxz[k] = -g_x_z_0_xxz[k] * ab_y + g_x_z_0_xxyz[k];

                g_x_z_y_xyy[k] = -g_x_z_0_xyy[k] * ab_y + g_x_z_0_xyyy[k];

                g_x_z_y_xyz[k] = -g_x_z_0_xyz[k] * ab_y + g_x_z_0_xyyz[k];

                g_x_z_y_xzz[k] = -g_x_z_0_xzz[k] * ab_y + g_x_z_0_xyzz[k];

                g_x_z_y_yyy[k] = -g_x_z_0_yyy[k] * ab_y + g_x_z_0_yyyy[k];

                g_x_z_y_yyz[k] = -g_x_z_0_yyz[k] * ab_y + g_x_z_0_yyyz[k];

                g_x_z_y_yzz[k] = -g_x_z_0_yzz[k] * ab_y + g_x_z_0_yyzz[k];

                g_x_z_y_zzz[k] = -g_x_z_0_zzz[k] * ab_y + g_x_z_0_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_x_z_z_xxx = cbuffer.data(pf_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_z_z_xxy = cbuffer.data(pf_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_z_z_xxz = cbuffer.data(pf_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_z_z_xyy = cbuffer.data(pf_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_z_z_xyz = cbuffer.data(pf_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_z_z_xzz = cbuffer.data(pf_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_z_z_yyy = cbuffer.data(pf_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_z_z_yyz = cbuffer.data(pf_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_z_z_yzz = cbuffer.data(pf_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_z_z_zzz = cbuffer.data(pf_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxx, g_x_0_0_xxy, g_x_0_0_xxz, g_x_0_0_xyy, g_x_0_0_xyz, g_x_0_0_xzz, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yzz, g_x_0_0_zzz, g_x_z_0_xxx, g_x_z_0_xxxz, g_x_z_0_xxy, g_x_z_0_xxyz, g_x_z_0_xxz, g_x_z_0_xxzz, g_x_z_0_xyy, g_x_z_0_xyyz, g_x_z_0_xyz, g_x_z_0_xyzz, g_x_z_0_xzz, g_x_z_0_xzzz, g_x_z_0_yyy, g_x_z_0_yyyz, g_x_z_0_yyz, g_x_z_0_yyzz, g_x_z_0_yzz, g_x_z_0_yzzz, g_x_z_0_zzz, g_x_z_0_zzzz, g_x_z_z_xxx, g_x_z_z_xxy, g_x_z_z_xxz, g_x_z_z_xyy, g_x_z_z_xyz, g_x_z_z_xzz, g_x_z_z_yyy, g_x_z_z_yyz, g_x_z_z_yzz, g_x_z_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_z_xxx[k] = g_x_0_0_xxx[k] - g_x_z_0_xxx[k] * ab_z + g_x_z_0_xxxz[k];

                g_x_z_z_xxy[k] = g_x_0_0_xxy[k] - g_x_z_0_xxy[k] * ab_z + g_x_z_0_xxyz[k];

                g_x_z_z_xxz[k] = g_x_0_0_xxz[k] - g_x_z_0_xxz[k] * ab_z + g_x_z_0_xxzz[k];

                g_x_z_z_xyy[k] = g_x_0_0_xyy[k] - g_x_z_0_xyy[k] * ab_z + g_x_z_0_xyyz[k];

                g_x_z_z_xyz[k] = g_x_0_0_xyz[k] - g_x_z_0_xyz[k] * ab_z + g_x_z_0_xyzz[k];

                g_x_z_z_xzz[k] = g_x_0_0_xzz[k] - g_x_z_0_xzz[k] * ab_z + g_x_z_0_xzzz[k];

                g_x_z_z_yyy[k] = g_x_0_0_yyy[k] - g_x_z_0_yyy[k] * ab_z + g_x_z_0_yyyz[k];

                g_x_z_z_yyz[k] = g_x_0_0_yyz[k] - g_x_z_0_yyz[k] * ab_z + g_x_z_0_yyzz[k];

                g_x_z_z_yzz[k] = g_x_0_0_yzz[k] - g_x_z_0_yzz[k] * ab_z + g_x_z_0_yzzz[k];

                g_x_z_z_zzz[k] = g_x_0_0_zzz[k] - g_x_z_0_zzz[k] * ab_z + g_x_z_0_zzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_y_x_x_xxx = cbuffer.data(pf_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_x_x_xxy = cbuffer.data(pf_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_x_x_xxz = cbuffer.data(pf_geom_11_off + 92 * ccomps * dcomps);

            auto g_y_x_x_xyy = cbuffer.data(pf_geom_11_off + 93 * ccomps * dcomps);

            auto g_y_x_x_xyz = cbuffer.data(pf_geom_11_off + 94 * ccomps * dcomps);

            auto g_y_x_x_xzz = cbuffer.data(pf_geom_11_off + 95 * ccomps * dcomps);

            auto g_y_x_x_yyy = cbuffer.data(pf_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_x_x_yyz = cbuffer.data(pf_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_x_x_yzz = cbuffer.data(pf_geom_11_off + 98 * ccomps * dcomps);

            auto g_y_x_x_zzz = cbuffer.data(pf_geom_11_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxx, g_y_0_0_xxy, g_y_0_0_xxz, g_y_0_0_xyy, g_y_0_0_xyz, g_y_0_0_xzz, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yzz, g_y_0_0_zzz, g_y_x_0_xxx, g_y_x_0_xxxx, g_y_x_0_xxxy, g_y_x_0_xxxz, g_y_x_0_xxy, g_y_x_0_xxyy, g_y_x_0_xxyz, g_y_x_0_xxz, g_y_x_0_xxzz, g_y_x_0_xyy, g_y_x_0_xyyy, g_y_x_0_xyyz, g_y_x_0_xyz, g_y_x_0_xyzz, g_y_x_0_xzz, g_y_x_0_xzzz, g_y_x_0_yyy, g_y_x_0_yyz, g_y_x_0_yzz, g_y_x_0_zzz, g_y_x_x_xxx, g_y_x_x_xxy, g_y_x_x_xxz, g_y_x_x_xyy, g_y_x_x_xyz, g_y_x_x_xzz, g_y_x_x_yyy, g_y_x_x_yyz, g_y_x_x_yzz, g_y_x_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_x_xxx[k] = g_y_0_0_xxx[k] - g_y_x_0_xxx[k] * ab_x + g_y_x_0_xxxx[k];

                g_y_x_x_xxy[k] = g_y_0_0_xxy[k] - g_y_x_0_xxy[k] * ab_x + g_y_x_0_xxxy[k];

                g_y_x_x_xxz[k] = g_y_0_0_xxz[k] - g_y_x_0_xxz[k] * ab_x + g_y_x_0_xxxz[k];

                g_y_x_x_xyy[k] = g_y_0_0_xyy[k] - g_y_x_0_xyy[k] * ab_x + g_y_x_0_xxyy[k];

                g_y_x_x_xyz[k] = g_y_0_0_xyz[k] - g_y_x_0_xyz[k] * ab_x + g_y_x_0_xxyz[k];

                g_y_x_x_xzz[k] = g_y_0_0_xzz[k] - g_y_x_0_xzz[k] * ab_x + g_y_x_0_xxzz[k];

                g_y_x_x_yyy[k] = g_y_0_0_yyy[k] - g_y_x_0_yyy[k] * ab_x + g_y_x_0_xyyy[k];

                g_y_x_x_yyz[k] = g_y_0_0_yyz[k] - g_y_x_0_yyz[k] * ab_x + g_y_x_0_xyyz[k];

                g_y_x_x_yzz[k] = g_y_0_0_yzz[k] - g_y_x_0_yzz[k] * ab_x + g_y_x_0_xyzz[k];

                g_y_x_x_zzz[k] = g_y_0_0_zzz[k] - g_y_x_0_zzz[k] * ab_x + g_y_x_0_xzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_y_x_y_xxx = cbuffer.data(pf_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_x_y_xxy = cbuffer.data(pf_geom_11_off + 101 * ccomps * dcomps);

            auto g_y_x_y_xxz = cbuffer.data(pf_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_x_y_xyy = cbuffer.data(pf_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_x_y_xyz = cbuffer.data(pf_geom_11_off + 104 * ccomps * dcomps);

            auto g_y_x_y_xzz = cbuffer.data(pf_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_x_y_yyy = cbuffer.data(pf_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_x_y_yyz = cbuffer.data(pf_geom_11_off + 107 * ccomps * dcomps);

            auto g_y_x_y_yzz = cbuffer.data(pf_geom_11_off + 108 * ccomps * dcomps);

            auto g_y_x_y_zzz = cbuffer.data(pf_geom_11_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxx, g_0_x_0_xxy, g_0_x_0_xxz, g_0_x_0_xyy, g_0_x_0_xyz, g_0_x_0_xzz, g_0_x_0_yyy, g_0_x_0_yyz, g_0_x_0_yzz, g_0_x_0_zzz, g_y_x_0_xxx, g_y_x_0_xxxy, g_y_x_0_xxy, g_y_x_0_xxyy, g_y_x_0_xxyz, g_y_x_0_xxz, g_y_x_0_xyy, g_y_x_0_xyyy, g_y_x_0_xyyz, g_y_x_0_xyz, g_y_x_0_xyzz, g_y_x_0_xzz, g_y_x_0_yyy, g_y_x_0_yyyy, g_y_x_0_yyyz, g_y_x_0_yyz, g_y_x_0_yyzz, g_y_x_0_yzz, g_y_x_0_yzzz, g_y_x_0_zzz, g_y_x_y_xxx, g_y_x_y_xxy, g_y_x_y_xxz, g_y_x_y_xyy, g_y_x_y_xyz, g_y_x_y_xzz, g_y_x_y_yyy, g_y_x_y_yyz, g_y_x_y_yzz, g_y_x_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_y_xxx[k] = -g_0_x_0_xxx[k] - g_y_x_0_xxx[k] * ab_y + g_y_x_0_xxxy[k];

                g_y_x_y_xxy[k] = -g_0_x_0_xxy[k] - g_y_x_0_xxy[k] * ab_y + g_y_x_0_xxyy[k];

                g_y_x_y_xxz[k] = -g_0_x_0_xxz[k] - g_y_x_0_xxz[k] * ab_y + g_y_x_0_xxyz[k];

                g_y_x_y_xyy[k] = -g_0_x_0_xyy[k] - g_y_x_0_xyy[k] * ab_y + g_y_x_0_xyyy[k];

                g_y_x_y_xyz[k] = -g_0_x_0_xyz[k] - g_y_x_0_xyz[k] * ab_y + g_y_x_0_xyyz[k];

                g_y_x_y_xzz[k] = -g_0_x_0_xzz[k] - g_y_x_0_xzz[k] * ab_y + g_y_x_0_xyzz[k];

                g_y_x_y_yyy[k] = -g_0_x_0_yyy[k] - g_y_x_0_yyy[k] * ab_y + g_y_x_0_yyyy[k];

                g_y_x_y_yyz[k] = -g_0_x_0_yyz[k] - g_y_x_0_yyz[k] * ab_y + g_y_x_0_yyyz[k];

                g_y_x_y_yzz[k] = -g_0_x_0_yzz[k] - g_y_x_0_yzz[k] * ab_y + g_y_x_0_yyzz[k];

                g_y_x_y_zzz[k] = -g_0_x_0_zzz[k] - g_y_x_0_zzz[k] * ab_y + g_y_x_0_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_y_x_z_xxx = cbuffer.data(pf_geom_11_off + 110 * ccomps * dcomps);

            auto g_y_x_z_xxy = cbuffer.data(pf_geom_11_off + 111 * ccomps * dcomps);

            auto g_y_x_z_xxz = cbuffer.data(pf_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_x_z_xyy = cbuffer.data(pf_geom_11_off + 113 * ccomps * dcomps);

            auto g_y_x_z_xyz = cbuffer.data(pf_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_x_z_xzz = cbuffer.data(pf_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_x_z_yyy = cbuffer.data(pf_geom_11_off + 116 * ccomps * dcomps);

            auto g_y_x_z_yyz = cbuffer.data(pf_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_x_z_yzz = cbuffer.data(pf_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_x_z_zzz = cbuffer.data(pf_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_0_xxx, g_y_x_0_xxxz, g_y_x_0_xxy, g_y_x_0_xxyz, g_y_x_0_xxz, g_y_x_0_xxzz, g_y_x_0_xyy, g_y_x_0_xyyz, g_y_x_0_xyz, g_y_x_0_xyzz, g_y_x_0_xzz, g_y_x_0_xzzz, g_y_x_0_yyy, g_y_x_0_yyyz, g_y_x_0_yyz, g_y_x_0_yyzz, g_y_x_0_yzz, g_y_x_0_yzzz, g_y_x_0_zzz, g_y_x_0_zzzz, g_y_x_z_xxx, g_y_x_z_xxy, g_y_x_z_xxz, g_y_x_z_xyy, g_y_x_z_xyz, g_y_x_z_xzz, g_y_x_z_yyy, g_y_x_z_yyz, g_y_x_z_yzz, g_y_x_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_z_xxx[k] = -g_y_x_0_xxx[k] * ab_z + g_y_x_0_xxxz[k];

                g_y_x_z_xxy[k] = -g_y_x_0_xxy[k] * ab_z + g_y_x_0_xxyz[k];

                g_y_x_z_xxz[k] = -g_y_x_0_xxz[k] * ab_z + g_y_x_0_xxzz[k];

                g_y_x_z_xyy[k] = -g_y_x_0_xyy[k] * ab_z + g_y_x_0_xyyz[k];

                g_y_x_z_xyz[k] = -g_y_x_0_xyz[k] * ab_z + g_y_x_0_xyzz[k];

                g_y_x_z_xzz[k] = -g_y_x_0_xzz[k] * ab_z + g_y_x_0_xzzz[k];

                g_y_x_z_yyy[k] = -g_y_x_0_yyy[k] * ab_z + g_y_x_0_yyyz[k];

                g_y_x_z_yyz[k] = -g_y_x_0_yyz[k] * ab_z + g_y_x_0_yyzz[k];

                g_y_x_z_yzz[k] = -g_y_x_0_yzz[k] * ab_z + g_y_x_0_yzzz[k];

                g_y_x_z_zzz[k] = -g_y_x_0_zzz[k] * ab_z + g_y_x_0_zzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_y_y_x_xxx = cbuffer.data(pf_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_y_x_xxy = cbuffer.data(pf_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_y_x_xxz = cbuffer.data(pf_geom_11_off + 122 * ccomps * dcomps);

            auto g_y_y_x_xyy = cbuffer.data(pf_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_y_x_xyz = cbuffer.data(pf_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_y_x_xzz = cbuffer.data(pf_geom_11_off + 125 * ccomps * dcomps);

            auto g_y_y_x_yyy = cbuffer.data(pf_geom_11_off + 126 * ccomps * dcomps);

            auto g_y_y_x_yyz = cbuffer.data(pf_geom_11_off + 127 * ccomps * dcomps);

            auto g_y_y_x_yzz = cbuffer.data(pf_geom_11_off + 128 * ccomps * dcomps);

            auto g_y_y_x_zzz = cbuffer.data(pf_geom_11_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_0_xxx, g_y_y_0_xxxx, g_y_y_0_xxxy, g_y_y_0_xxxz, g_y_y_0_xxy, g_y_y_0_xxyy, g_y_y_0_xxyz, g_y_y_0_xxz, g_y_y_0_xxzz, g_y_y_0_xyy, g_y_y_0_xyyy, g_y_y_0_xyyz, g_y_y_0_xyz, g_y_y_0_xyzz, g_y_y_0_xzz, g_y_y_0_xzzz, g_y_y_0_yyy, g_y_y_0_yyz, g_y_y_0_yzz, g_y_y_0_zzz, g_y_y_x_xxx, g_y_y_x_xxy, g_y_y_x_xxz, g_y_y_x_xyy, g_y_y_x_xyz, g_y_y_x_xzz, g_y_y_x_yyy, g_y_y_x_yyz, g_y_y_x_yzz, g_y_y_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_x_xxx[k] = -g_y_y_0_xxx[k] * ab_x + g_y_y_0_xxxx[k];

                g_y_y_x_xxy[k] = -g_y_y_0_xxy[k] * ab_x + g_y_y_0_xxxy[k];

                g_y_y_x_xxz[k] = -g_y_y_0_xxz[k] * ab_x + g_y_y_0_xxxz[k];

                g_y_y_x_xyy[k] = -g_y_y_0_xyy[k] * ab_x + g_y_y_0_xxyy[k];

                g_y_y_x_xyz[k] = -g_y_y_0_xyz[k] * ab_x + g_y_y_0_xxyz[k];

                g_y_y_x_xzz[k] = -g_y_y_0_xzz[k] * ab_x + g_y_y_0_xxzz[k];

                g_y_y_x_yyy[k] = -g_y_y_0_yyy[k] * ab_x + g_y_y_0_xyyy[k];

                g_y_y_x_yyz[k] = -g_y_y_0_yyz[k] * ab_x + g_y_y_0_xyyz[k];

                g_y_y_x_yzz[k] = -g_y_y_0_yzz[k] * ab_x + g_y_y_0_xyzz[k];

                g_y_y_x_zzz[k] = -g_y_y_0_zzz[k] * ab_x + g_y_y_0_xzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_y_y_y_xxx = cbuffer.data(pf_geom_11_off + 130 * ccomps * dcomps);

            auto g_y_y_y_xxy = cbuffer.data(pf_geom_11_off + 131 * ccomps * dcomps);

            auto g_y_y_y_xxz = cbuffer.data(pf_geom_11_off + 132 * ccomps * dcomps);

            auto g_y_y_y_xyy = cbuffer.data(pf_geom_11_off + 133 * ccomps * dcomps);

            auto g_y_y_y_xyz = cbuffer.data(pf_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_y_y_xzz = cbuffer.data(pf_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_y_y_yyy = cbuffer.data(pf_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_y_y_yyz = cbuffer.data(pf_geom_11_off + 137 * ccomps * dcomps);

            auto g_y_y_y_yzz = cbuffer.data(pf_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_y_y_zzz = cbuffer.data(pf_geom_11_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxx, g_0_y_0_xxy, g_0_y_0_xxz, g_0_y_0_xyy, g_0_y_0_xyz, g_0_y_0_xzz, g_0_y_0_yyy, g_0_y_0_yyz, g_0_y_0_yzz, g_0_y_0_zzz, g_y_0_0_xxx, g_y_0_0_xxy, g_y_0_0_xxz, g_y_0_0_xyy, g_y_0_0_xyz, g_y_0_0_xzz, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yzz, g_y_0_0_zzz, g_y_y_0_xxx, g_y_y_0_xxxy, g_y_y_0_xxy, g_y_y_0_xxyy, g_y_y_0_xxyz, g_y_y_0_xxz, g_y_y_0_xyy, g_y_y_0_xyyy, g_y_y_0_xyyz, g_y_y_0_xyz, g_y_y_0_xyzz, g_y_y_0_xzz, g_y_y_0_yyy, g_y_y_0_yyyy, g_y_y_0_yyyz, g_y_y_0_yyz, g_y_y_0_yyzz, g_y_y_0_yzz, g_y_y_0_yzzz, g_y_y_0_zzz, g_y_y_y_xxx, g_y_y_y_xxy, g_y_y_y_xxz, g_y_y_y_xyy, g_y_y_y_xyz, g_y_y_y_xzz, g_y_y_y_yyy, g_y_y_y_yyz, g_y_y_y_yzz, g_y_y_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_y_xxx[k] = -g_0_y_0_xxx[k] + g_y_0_0_xxx[k] - g_y_y_0_xxx[k] * ab_y + g_y_y_0_xxxy[k];

                g_y_y_y_xxy[k] = -g_0_y_0_xxy[k] + g_y_0_0_xxy[k] - g_y_y_0_xxy[k] * ab_y + g_y_y_0_xxyy[k];

                g_y_y_y_xxz[k] = -g_0_y_0_xxz[k] + g_y_0_0_xxz[k] - g_y_y_0_xxz[k] * ab_y + g_y_y_0_xxyz[k];

                g_y_y_y_xyy[k] = -g_0_y_0_xyy[k] + g_y_0_0_xyy[k] - g_y_y_0_xyy[k] * ab_y + g_y_y_0_xyyy[k];

                g_y_y_y_xyz[k] = -g_0_y_0_xyz[k] + g_y_0_0_xyz[k] - g_y_y_0_xyz[k] * ab_y + g_y_y_0_xyyz[k];

                g_y_y_y_xzz[k] = -g_0_y_0_xzz[k] + g_y_0_0_xzz[k] - g_y_y_0_xzz[k] * ab_y + g_y_y_0_xyzz[k];

                g_y_y_y_yyy[k] = -g_0_y_0_yyy[k] + g_y_0_0_yyy[k] - g_y_y_0_yyy[k] * ab_y + g_y_y_0_yyyy[k];

                g_y_y_y_yyz[k] = -g_0_y_0_yyz[k] + g_y_0_0_yyz[k] - g_y_y_0_yyz[k] * ab_y + g_y_y_0_yyyz[k];

                g_y_y_y_yzz[k] = -g_0_y_0_yzz[k] + g_y_0_0_yzz[k] - g_y_y_0_yzz[k] * ab_y + g_y_y_0_yyzz[k];

                g_y_y_y_zzz[k] = -g_0_y_0_zzz[k] + g_y_0_0_zzz[k] - g_y_y_0_zzz[k] * ab_y + g_y_y_0_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_y_y_z_xxx = cbuffer.data(pf_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_y_z_xxy = cbuffer.data(pf_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_y_z_xxz = cbuffer.data(pf_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_y_z_xyy = cbuffer.data(pf_geom_11_off + 143 * ccomps * dcomps);

            auto g_y_y_z_xyz = cbuffer.data(pf_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_y_z_xzz = cbuffer.data(pf_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_y_z_yyy = cbuffer.data(pf_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_y_z_yyz = cbuffer.data(pf_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_y_z_yzz = cbuffer.data(pf_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_y_z_zzz = cbuffer.data(pf_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_0_xxx, g_y_y_0_xxxz, g_y_y_0_xxy, g_y_y_0_xxyz, g_y_y_0_xxz, g_y_y_0_xxzz, g_y_y_0_xyy, g_y_y_0_xyyz, g_y_y_0_xyz, g_y_y_0_xyzz, g_y_y_0_xzz, g_y_y_0_xzzz, g_y_y_0_yyy, g_y_y_0_yyyz, g_y_y_0_yyz, g_y_y_0_yyzz, g_y_y_0_yzz, g_y_y_0_yzzz, g_y_y_0_zzz, g_y_y_0_zzzz, g_y_y_z_xxx, g_y_y_z_xxy, g_y_y_z_xxz, g_y_y_z_xyy, g_y_y_z_xyz, g_y_y_z_xzz, g_y_y_z_yyy, g_y_y_z_yyz, g_y_y_z_yzz, g_y_y_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_z_xxx[k] = -g_y_y_0_xxx[k] * ab_z + g_y_y_0_xxxz[k];

                g_y_y_z_xxy[k] = -g_y_y_0_xxy[k] * ab_z + g_y_y_0_xxyz[k];

                g_y_y_z_xxz[k] = -g_y_y_0_xxz[k] * ab_z + g_y_y_0_xxzz[k];

                g_y_y_z_xyy[k] = -g_y_y_0_xyy[k] * ab_z + g_y_y_0_xyyz[k];

                g_y_y_z_xyz[k] = -g_y_y_0_xyz[k] * ab_z + g_y_y_0_xyzz[k];

                g_y_y_z_xzz[k] = -g_y_y_0_xzz[k] * ab_z + g_y_y_0_xzzz[k];

                g_y_y_z_yyy[k] = -g_y_y_0_yyy[k] * ab_z + g_y_y_0_yyyz[k];

                g_y_y_z_yyz[k] = -g_y_y_0_yyz[k] * ab_z + g_y_y_0_yyzz[k];

                g_y_y_z_yzz[k] = -g_y_y_0_yzz[k] * ab_z + g_y_y_0_yzzz[k];

                g_y_y_z_zzz[k] = -g_y_y_0_zzz[k] * ab_z + g_y_y_0_zzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_y_z_x_xxx = cbuffer.data(pf_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_z_x_xxy = cbuffer.data(pf_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_z_x_xxz = cbuffer.data(pf_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_z_x_xyy = cbuffer.data(pf_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_z_x_xyz = cbuffer.data(pf_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_z_x_xzz = cbuffer.data(pf_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_z_x_yyy = cbuffer.data(pf_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_z_x_yyz = cbuffer.data(pf_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_z_x_yzz = cbuffer.data(pf_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_z_x_zzz = cbuffer.data(pf_geom_11_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_0_xxx, g_y_z_0_xxxx, g_y_z_0_xxxy, g_y_z_0_xxxz, g_y_z_0_xxy, g_y_z_0_xxyy, g_y_z_0_xxyz, g_y_z_0_xxz, g_y_z_0_xxzz, g_y_z_0_xyy, g_y_z_0_xyyy, g_y_z_0_xyyz, g_y_z_0_xyz, g_y_z_0_xyzz, g_y_z_0_xzz, g_y_z_0_xzzz, g_y_z_0_yyy, g_y_z_0_yyz, g_y_z_0_yzz, g_y_z_0_zzz, g_y_z_x_xxx, g_y_z_x_xxy, g_y_z_x_xxz, g_y_z_x_xyy, g_y_z_x_xyz, g_y_z_x_xzz, g_y_z_x_yyy, g_y_z_x_yyz, g_y_z_x_yzz, g_y_z_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_x_xxx[k] = -g_y_z_0_xxx[k] * ab_x + g_y_z_0_xxxx[k];

                g_y_z_x_xxy[k] = -g_y_z_0_xxy[k] * ab_x + g_y_z_0_xxxy[k];

                g_y_z_x_xxz[k] = -g_y_z_0_xxz[k] * ab_x + g_y_z_0_xxxz[k];

                g_y_z_x_xyy[k] = -g_y_z_0_xyy[k] * ab_x + g_y_z_0_xxyy[k];

                g_y_z_x_xyz[k] = -g_y_z_0_xyz[k] * ab_x + g_y_z_0_xxyz[k];

                g_y_z_x_xzz[k] = -g_y_z_0_xzz[k] * ab_x + g_y_z_0_xxzz[k];

                g_y_z_x_yyy[k] = -g_y_z_0_yyy[k] * ab_x + g_y_z_0_xyyy[k];

                g_y_z_x_yyz[k] = -g_y_z_0_yyz[k] * ab_x + g_y_z_0_xyyz[k];

                g_y_z_x_yzz[k] = -g_y_z_0_yzz[k] * ab_x + g_y_z_0_xyzz[k];

                g_y_z_x_zzz[k] = -g_y_z_0_zzz[k] * ab_x + g_y_z_0_xzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_y_z_y_xxx = cbuffer.data(pf_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_z_y_xxy = cbuffer.data(pf_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_z_y_xxz = cbuffer.data(pf_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_z_y_xyy = cbuffer.data(pf_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_z_y_xyz = cbuffer.data(pf_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_z_y_xzz = cbuffer.data(pf_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_z_y_yyy = cbuffer.data(pf_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_z_y_yyz = cbuffer.data(pf_geom_11_off + 167 * ccomps * dcomps);

            auto g_y_z_y_yzz = cbuffer.data(pf_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_z_y_zzz = cbuffer.data(pf_geom_11_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxx, g_0_z_0_xxy, g_0_z_0_xxz, g_0_z_0_xyy, g_0_z_0_xyz, g_0_z_0_xzz, g_0_z_0_yyy, g_0_z_0_yyz, g_0_z_0_yzz, g_0_z_0_zzz, g_y_z_0_xxx, g_y_z_0_xxxy, g_y_z_0_xxy, g_y_z_0_xxyy, g_y_z_0_xxyz, g_y_z_0_xxz, g_y_z_0_xyy, g_y_z_0_xyyy, g_y_z_0_xyyz, g_y_z_0_xyz, g_y_z_0_xyzz, g_y_z_0_xzz, g_y_z_0_yyy, g_y_z_0_yyyy, g_y_z_0_yyyz, g_y_z_0_yyz, g_y_z_0_yyzz, g_y_z_0_yzz, g_y_z_0_yzzz, g_y_z_0_zzz, g_y_z_y_xxx, g_y_z_y_xxy, g_y_z_y_xxz, g_y_z_y_xyy, g_y_z_y_xyz, g_y_z_y_xzz, g_y_z_y_yyy, g_y_z_y_yyz, g_y_z_y_yzz, g_y_z_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_y_xxx[k] = -g_0_z_0_xxx[k] - g_y_z_0_xxx[k] * ab_y + g_y_z_0_xxxy[k];

                g_y_z_y_xxy[k] = -g_0_z_0_xxy[k] - g_y_z_0_xxy[k] * ab_y + g_y_z_0_xxyy[k];

                g_y_z_y_xxz[k] = -g_0_z_0_xxz[k] - g_y_z_0_xxz[k] * ab_y + g_y_z_0_xxyz[k];

                g_y_z_y_xyy[k] = -g_0_z_0_xyy[k] - g_y_z_0_xyy[k] * ab_y + g_y_z_0_xyyy[k];

                g_y_z_y_xyz[k] = -g_0_z_0_xyz[k] - g_y_z_0_xyz[k] * ab_y + g_y_z_0_xyyz[k];

                g_y_z_y_xzz[k] = -g_0_z_0_xzz[k] - g_y_z_0_xzz[k] * ab_y + g_y_z_0_xyzz[k];

                g_y_z_y_yyy[k] = -g_0_z_0_yyy[k] - g_y_z_0_yyy[k] * ab_y + g_y_z_0_yyyy[k];

                g_y_z_y_yyz[k] = -g_0_z_0_yyz[k] - g_y_z_0_yyz[k] * ab_y + g_y_z_0_yyyz[k];

                g_y_z_y_yzz[k] = -g_0_z_0_yzz[k] - g_y_z_0_yzz[k] * ab_y + g_y_z_0_yyzz[k];

                g_y_z_y_zzz[k] = -g_0_z_0_zzz[k] - g_y_z_0_zzz[k] * ab_y + g_y_z_0_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_y_z_z_xxx = cbuffer.data(pf_geom_11_off + 170 * ccomps * dcomps);

            auto g_y_z_z_xxy = cbuffer.data(pf_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_z_z_xxz = cbuffer.data(pf_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_z_z_xyy = cbuffer.data(pf_geom_11_off + 173 * ccomps * dcomps);

            auto g_y_z_z_xyz = cbuffer.data(pf_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_z_z_xzz = cbuffer.data(pf_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_z_z_yyy = cbuffer.data(pf_geom_11_off + 176 * ccomps * dcomps);

            auto g_y_z_z_yyz = cbuffer.data(pf_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_z_z_yzz = cbuffer.data(pf_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_z_z_zzz = cbuffer.data(pf_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxx, g_y_0_0_xxy, g_y_0_0_xxz, g_y_0_0_xyy, g_y_0_0_xyz, g_y_0_0_xzz, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yzz, g_y_0_0_zzz, g_y_z_0_xxx, g_y_z_0_xxxz, g_y_z_0_xxy, g_y_z_0_xxyz, g_y_z_0_xxz, g_y_z_0_xxzz, g_y_z_0_xyy, g_y_z_0_xyyz, g_y_z_0_xyz, g_y_z_0_xyzz, g_y_z_0_xzz, g_y_z_0_xzzz, g_y_z_0_yyy, g_y_z_0_yyyz, g_y_z_0_yyz, g_y_z_0_yyzz, g_y_z_0_yzz, g_y_z_0_yzzz, g_y_z_0_zzz, g_y_z_0_zzzz, g_y_z_z_xxx, g_y_z_z_xxy, g_y_z_z_xxz, g_y_z_z_xyy, g_y_z_z_xyz, g_y_z_z_xzz, g_y_z_z_yyy, g_y_z_z_yyz, g_y_z_z_yzz, g_y_z_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_z_xxx[k] = g_y_0_0_xxx[k] - g_y_z_0_xxx[k] * ab_z + g_y_z_0_xxxz[k];

                g_y_z_z_xxy[k] = g_y_0_0_xxy[k] - g_y_z_0_xxy[k] * ab_z + g_y_z_0_xxyz[k];

                g_y_z_z_xxz[k] = g_y_0_0_xxz[k] - g_y_z_0_xxz[k] * ab_z + g_y_z_0_xxzz[k];

                g_y_z_z_xyy[k] = g_y_0_0_xyy[k] - g_y_z_0_xyy[k] * ab_z + g_y_z_0_xyyz[k];

                g_y_z_z_xyz[k] = g_y_0_0_xyz[k] - g_y_z_0_xyz[k] * ab_z + g_y_z_0_xyzz[k];

                g_y_z_z_xzz[k] = g_y_0_0_xzz[k] - g_y_z_0_xzz[k] * ab_z + g_y_z_0_xzzz[k];

                g_y_z_z_yyy[k] = g_y_0_0_yyy[k] - g_y_z_0_yyy[k] * ab_z + g_y_z_0_yyyz[k];

                g_y_z_z_yyz[k] = g_y_0_0_yyz[k] - g_y_z_0_yyz[k] * ab_z + g_y_z_0_yyzz[k];

                g_y_z_z_yzz[k] = g_y_0_0_yzz[k] - g_y_z_0_yzz[k] * ab_z + g_y_z_0_yzzz[k];

                g_y_z_z_zzz[k] = g_y_0_0_zzz[k] - g_y_z_0_zzz[k] * ab_z + g_y_z_0_zzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_z_x_x_xxx = cbuffer.data(pf_geom_11_off + 180 * ccomps * dcomps);

            auto g_z_x_x_xxy = cbuffer.data(pf_geom_11_off + 181 * ccomps * dcomps);

            auto g_z_x_x_xxz = cbuffer.data(pf_geom_11_off + 182 * ccomps * dcomps);

            auto g_z_x_x_xyy = cbuffer.data(pf_geom_11_off + 183 * ccomps * dcomps);

            auto g_z_x_x_xyz = cbuffer.data(pf_geom_11_off + 184 * ccomps * dcomps);

            auto g_z_x_x_xzz = cbuffer.data(pf_geom_11_off + 185 * ccomps * dcomps);

            auto g_z_x_x_yyy = cbuffer.data(pf_geom_11_off + 186 * ccomps * dcomps);

            auto g_z_x_x_yyz = cbuffer.data(pf_geom_11_off + 187 * ccomps * dcomps);

            auto g_z_x_x_yzz = cbuffer.data(pf_geom_11_off + 188 * ccomps * dcomps);

            auto g_z_x_x_zzz = cbuffer.data(pf_geom_11_off + 189 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxx, g_z_0_0_xxy, g_z_0_0_xxz, g_z_0_0_xyy, g_z_0_0_xyz, g_z_0_0_xzz, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yzz, g_z_0_0_zzz, g_z_x_0_xxx, g_z_x_0_xxxx, g_z_x_0_xxxy, g_z_x_0_xxxz, g_z_x_0_xxy, g_z_x_0_xxyy, g_z_x_0_xxyz, g_z_x_0_xxz, g_z_x_0_xxzz, g_z_x_0_xyy, g_z_x_0_xyyy, g_z_x_0_xyyz, g_z_x_0_xyz, g_z_x_0_xyzz, g_z_x_0_xzz, g_z_x_0_xzzz, g_z_x_0_yyy, g_z_x_0_yyz, g_z_x_0_yzz, g_z_x_0_zzz, g_z_x_x_xxx, g_z_x_x_xxy, g_z_x_x_xxz, g_z_x_x_xyy, g_z_x_x_xyz, g_z_x_x_xzz, g_z_x_x_yyy, g_z_x_x_yyz, g_z_x_x_yzz, g_z_x_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_x_xxx[k] = g_z_0_0_xxx[k] - g_z_x_0_xxx[k] * ab_x + g_z_x_0_xxxx[k];

                g_z_x_x_xxy[k] = g_z_0_0_xxy[k] - g_z_x_0_xxy[k] * ab_x + g_z_x_0_xxxy[k];

                g_z_x_x_xxz[k] = g_z_0_0_xxz[k] - g_z_x_0_xxz[k] * ab_x + g_z_x_0_xxxz[k];

                g_z_x_x_xyy[k] = g_z_0_0_xyy[k] - g_z_x_0_xyy[k] * ab_x + g_z_x_0_xxyy[k];

                g_z_x_x_xyz[k] = g_z_0_0_xyz[k] - g_z_x_0_xyz[k] * ab_x + g_z_x_0_xxyz[k];

                g_z_x_x_xzz[k] = g_z_0_0_xzz[k] - g_z_x_0_xzz[k] * ab_x + g_z_x_0_xxzz[k];

                g_z_x_x_yyy[k] = g_z_0_0_yyy[k] - g_z_x_0_yyy[k] * ab_x + g_z_x_0_xyyy[k];

                g_z_x_x_yyz[k] = g_z_0_0_yyz[k] - g_z_x_0_yyz[k] * ab_x + g_z_x_0_xyyz[k];

                g_z_x_x_yzz[k] = g_z_0_0_yzz[k] - g_z_x_0_yzz[k] * ab_x + g_z_x_0_xyzz[k];

                g_z_x_x_zzz[k] = g_z_0_0_zzz[k] - g_z_x_0_zzz[k] * ab_x + g_z_x_0_xzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_z_x_y_xxx = cbuffer.data(pf_geom_11_off + 190 * ccomps * dcomps);

            auto g_z_x_y_xxy = cbuffer.data(pf_geom_11_off + 191 * ccomps * dcomps);

            auto g_z_x_y_xxz = cbuffer.data(pf_geom_11_off + 192 * ccomps * dcomps);

            auto g_z_x_y_xyy = cbuffer.data(pf_geom_11_off + 193 * ccomps * dcomps);

            auto g_z_x_y_xyz = cbuffer.data(pf_geom_11_off + 194 * ccomps * dcomps);

            auto g_z_x_y_xzz = cbuffer.data(pf_geom_11_off + 195 * ccomps * dcomps);

            auto g_z_x_y_yyy = cbuffer.data(pf_geom_11_off + 196 * ccomps * dcomps);

            auto g_z_x_y_yyz = cbuffer.data(pf_geom_11_off + 197 * ccomps * dcomps);

            auto g_z_x_y_yzz = cbuffer.data(pf_geom_11_off + 198 * ccomps * dcomps);

            auto g_z_x_y_zzz = cbuffer.data(pf_geom_11_off + 199 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_0_xxx, g_z_x_0_xxxy, g_z_x_0_xxy, g_z_x_0_xxyy, g_z_x_0_xxyz, g_z_x_0_xxz, g_z_x_0_xyy, g_z_x_0_xyyy, g_z_x_0_xyyz, g_z_x_0_xyz, g_z_x_0_xyzz, g_z_x_0_xzz, g_z_x_0_yyy, g_z_x_0_yyyy, g_z_x_0_yyyz, g_z_x_0_yyz, g_z_x_0_yyzz, g_z_x_0_yzz, g_z_x_0_yzzz, g_z_x_0_zzz, g_z_x_y_xxx, g_z_x_y_xxy, g_z_x_y_xxz, g_z_x_y_xyy, g_z_x_y_xyz, g_z_x_y_xzz, g_z_x_y_yyy, g_z_x_y_yyz, g_z_x_y_yzz, g_z_x_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_y_xxx[k] = -g_z_x_0_xxx[k] * ab_y + g_z_x_0_xxxy[k];

                g_z_x_y_xxy[k] = -g_z_x_0_xxy[k] * ab_y + g_z_x_0_xxyy[k];

                g_z_x_y_xxz[k] = -g_z_x_0_xxz[k] * ab_y + g_z_x_0_xxyz[k];

                g_z_x_y_xyy[k] = -g_z_x_0_xyy[k] * ab_y + g_z_x_0_xyyy[k];

                g_z_x_y_xyz[k] = -g_z_x_0_xyz[k] * ab_y + g_z_x_0_xyyz[k];

                g_z_x_y_xzz[k] = -g_z_x_0_xzz[k] * ab_y + g_z_x_0_xyzz[k];

                g_z_x_y_yyy[k] = -g_z_x_0_yyy[k] * ab_y + g_z_x_0_yyyy[k];

                g_z_x_y_yyz[k] = -g_z_x_0_yyz[k] * ab_y + g_z_x_0_yyyz[k];

                g_z_x_y_yzz[k] = -g_z_x_0_yzz[k] * ab_y + g_z_x_0_yyzz[k];

                g_z_x_y_zzz[k] = -g_z_x_0_zzz[k] * ab_y + g_z_x_0_yzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_z_x_z_xxx = cbuffer.data(pf_geom_11_off + 200 * ccomps * dcomps);

            auto g_z_x_z_xxy = cbuffer.data(pf_geom_11_off + 201 * ccomps * dcomps);

            auto g_z_x_z_xxz = cbuffer.data(pf_geom_11_off + 202 * ccomps * dcomps);

            auto g_z_x_z_xyy = cbuffer.data(pf_geom_11_off + 203 * ccomps * dcomps);

            auto g_z_x_z_xyz = cbuffer.data(pf_geom_11_off + 204 * ccomps * dcomps);

            auto g_z_x_z_xzz = cbuffer.data(pf_geom_11_off + 205 * ccomps * dcomps);

            auto g_z_x_z_yyy = cbuffer.data(pf_geom_11_off + 206 * ccomps * dcomps);

            auto g_z_x_z_yyz = cbuffer.data(pf_geom_11_off + 207 * ccomps * dcomps);

            auto g_z_x_z_yzz = cbuffer.data(pf_geom_11_off + 208 * ccomps * dcomps);

            auto g_z_x_z_zzz = cbuffer.data(pf_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxx, g_0_x_0_xxy, g_0_x_0_xxz, g_0_x_0_xyy, g_0_x_0_xyz, g_0_x_0_xzz, g_0_x_0_yyy, g_0_x_0_yyz, g_0_x_0_yzz, g_0_x_0_zzz, g_z_x_0_xxx, g_z_x_0_xxxz, g_z_x_0_xxy, g_z_x_0_xxyz, g_z_x_0_xxz, g_z_x_0_xxzz, g_z_x_0_xyy, g_z_x_0_xyyz, g_z_x_0_xyz, g_z_x_0_xyzz, g_z_x_0_xzz, g_z_x_0_xzzz, g_z_x_0_yyy, g_z_x_0_yyyz, g_z_x_0_yyz, g_z_x_0_yyzz, g_z_x_0_yzz, g_z_x_0_yzzz, g_z_x_0_zzz, g_z_x_0_zzzz, g_z_x_z_xxx, g_z_x_z_xxy, g_z_x_z_xxz, g_z_x_z_xyy, g_z_x_z_xyz, g_z_x_z_xzz, g_z_x_z_yyy, g_z_x_z_yyz, g_z_x_z_yzz, g_z_x_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_z_xxx[k] = -g_0_x_0_xxx[k] - g_z_x_0_xxx[k] * ab_z + g_z_x_0_xxxz[k];

                g_z_x_z_xxy[k] = -g_0_x_0_xxy[k] - g_z_x_0_xxy[k] * ab_z + g_z_x_0_xxyz[k];

                g_z_x_z_xxz[k] = -g_0_x_0_xxz[k] - g_z_x_0_xxz[k] * ab_z + g_z_x_0_xxzz[k];

                g_z_x_z_xyy[k] = -g_0_x_0_xyy[k] - g_z_x_0_xyy[k] * ab_z + g_z_x_0_xyyz[k];

                g_z_x_z_xyz[k] = -g_0_x_0_xyz[k] - g_z_x_0_xyz[k] * ab_z + g_z_x_0_xyzz[k];

                g_z_x_z_xzz[k] = -g_0_x_0_xzz[k] - g_z_x_0_xzz[k] * ab_z + g_z_x_0_xzzz[k];

                g_z_x_z_yyy[k] = -g_0_x_0_yyy[k] - g_z_x_0_yyy[k] * ab_z + g_z_x_0_yyyz[k];

                g_z_x_z_yyz[k] = -g_0_x_0_yyz[k] - g_z_x_0_yyz[k] * ab_z + g_z_x_0_yyzz[k];

                g_z_x_z_yzz[k] = -g_0_x_0_yzz[k] - g_z_x_0_yzz[k] * ab_z + g_z_x_0_yzzz[k];

                g_z_x_z_zzz[k] = -g_0_x_0_zzz[k] - g_z_x_0_zzz[k] * ab_z + g_z_x_0_zzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_z_y_x_xxx = cbuffer.data(pf_geom_11_off + 210 * ccomps * dcomps);

            auto g_z_y_x_xxy = cbuffer.data(pf_geom_11_off + 211 * ccomps * dcomps);

            auto g_z_y_x_xxz = cbuffer.data(pf_geom_11_off + 212 * ccomps * dcomps);

            auto g_z_y_x_xyy = cbuffer.data(pf_geom_11_off + 213 * ccomps * dcomps);

            auto g_z_y_x_xyz = cbuffer.data(pf_geom_11_off + 214 * ccomps * dcomps);

            auto g_z_y_x_xzz = cbuffer.data(pf_geom_11_off + 215 * ccomps * dcomps);

            auto g_z_y_x_yyy = cbuffer.data(pf_geom_11_off + 216 * ccomps * dcomps);

            auto g_z_y_x_yyz = cbuffer.data(pf_geom_11_off + 217 * ccomps * dcomps);

            auto g_z_y_x_yzz = cbuffer.data(pf_geom_11_off + 218 * ccomps * dcomps);

            auto g_z_y_x_zzz = cbuffer.data(pf_geom_11_off + 219 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_0_xxx, g_z_y_0_xxxx, g_z_y_0_xxxy, g_z_y_0_xxxz, g_z_y_0_xxy, g_z_y_0_xxyy, g_z_y_0_xxyz, g_z_y_0_xxz, g_z_y_0_xxzz, g_z_y_0_xyy, g_z_y_0_xyyy, g_z_y_0_xyyz, g_z_y_0_xyz, g_z_y_0_xyzz, g_z_y_0_xzz, g_z_y_0_xzzz, g_z_y_0_yyy, g_z_y_0_yyz, g_z_y_0_yzz, g_z_y_0_zzz, g_z_y_x_xxx, g_z_y_x_xxy, g_z_y_x_xxz, g_z_y_x_xyy, g_z_y_x_xyz, g_z_y_x_xzz, g_z_y_x_yyy, g_z_y_x_yyz, g_z_y_x_yzz, g_z_y_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_x_xxx[k] = -g_z_y_0_xxx[k] * ab_x + g_z_y_0_xxxx[k];

                g_z_y_x_xxy[k] = -g_z_y_0_xxy[k] * ab_x + g_z_y_0_xxxy[k];

                g_z_y_x_xxz[k] = -g_z_y_0_xxz[k] * ab_x + g_z_y_0_xxxz[k];

                g_z_y_x_xyy[k] = -g_z_y_0_xyy[k] * ab_x + g_z_y_0_xxyy[k];

                g_z_y_x_xyz[k] = -g_z_y_0_xyz[k] * ab_x + g_z_y_0_xxyz[k];

                g_z_y_x_xzz[k] = -g_z_y_0_xzz[k] * ab_x + g_z_y_0_xxzz[k];

                g_z_y_x_yyy[k] = -g_z_y_0_yyy[k] * ab_x + g_z_y_0_xyyy[k];

                g_z_y_x_yyz[k] = -g_z_y_0_yyz[k] * ab_x + g_z_y_0_xyyz[k];

                g_z_y_x_yzz[k] = -g_z_y_0_yzz[k] * ab_x + g_z_y_0_xyzz[k];

                g_z_y_x_zzz[k] = -g_z_y_0_zzz[k] * ab_x + g_z_y_0_xzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_z_y_y_xxx = cbuffer.data(pf_geom_11_off + 220 * ccomps * dcomps);

            auto g_z_y_y_xxy = cbuffer.data(pf_geom_11_off + 221 * ccomps * dcomps);

            auto g_z_y_y_xxz = cbuffer.data(pf_geom_11_off + 222 * ccomps * dcomps);

            auto g_z_y_y_xyy = cbuffer.data(pf_geom_11_off + 223 * ccomps * dcomps);

            auto g_z_y_y_xyz = cbuffer.data(pf_geom_11_off + 224 * ccomps * dcomps);

            auto g_z_y_y_xzz = cbuffer.data(pf_geom_11_off + 225 * ccomps * dcomps);

            auto g_z_y_y_yyy = cbuffer.data(pf_geom_11_off + 226 * ccomps * dcomps);

            auto g_z_y_y_yyz = cbuffer.data(pf_geom_11_off + 227 * ccomps * dcomps);

            auto g_z_y_y_yzz = cbuffer.data(pf_geom_11_off + 228 * ccomps * dcomps);

            auto g_z_y_y_zzz = cbuffer.data(pf_geom_11_off + 229 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxx, g_z_0_0_xxy, g_z_0_0_xxz, g_z_0_0_xyy, g_z_0_0_xyz, g_z_0_0_xzz, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yzz, g_z_0_0_zzz, g_z_y_0_xxx, g_z_y_0_xxxy, g_z_y_0_xxy, g_z_y_0_xxyy, g_z_y_0_xxyz, g_z_y_0_xxz, g_z_y_0_xyy, g_z_y_0_xyyy, g_z_y_0_xyyz, g_z_y_0_xyz, g_z_y_0_xyzz, g_z_y_0_xzz, g_z_y_0_yyy, g_z_y_0_yyyy, g_z_y_0_yyyz, g_z_y_0_yyz, g_z_y_0_yyzz, g_z_y_0_yzz, g_z_y_0_yzzz, g_z_y_0_zzz, g_z_y_y_xxx, g_z_y_y_xxy, g_z_y_y_xxz, g_z_y_y_xyy, g_z_y_y_xyz, g_z_y_y_xzz, g_z_y_y_yyy, g_z_y_y_yyz, g_z_y_y_yzz, g_z_y_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_y_xxx[k] = g_z_0_0_xxx[k] - g_z_y_0_xxx[k] * ab_y + g_z_y_0_xxxy[k];

                g_z_y_y_xxy[k] = g_z_0_0_xxy[k] - g_z_y_0_xxy[k] * ab_y + g_z_y_0_xxyy[k];

                g_z_y_y_xxz[k] = g_z_0_0_xxz[k] - g_z_y_0_xxz[k] * ab_y + g_z_y_0_xxyz[k];

                g_z_y_y_xyy[k] = g_z_0_0_xyy[k] - g_z_y_0_xyy[k] * ab_y + g_z_y_0_xyyy[k];

                g_z_y_y_xyz[k] = g_z_0_0_xyz[k] - g_z_y_0_xyz[k] * ab_y + g_z_y_0_xyyz[k];

                g_z_y_y_xzz[k] = g_z_0_0_xzz[k] - g_z_y_0_xzz[k] * ab_y + g_z_y_0_xyzz[k];

                g_z_y_y_yyy[k] = g_z_0_0_yyy[k] - g_z_y_0_yyy[k] * ab_y + g_z_y_0_yyyy[k];

                g_z_y_y_yyz[k] = g_z_0_0_yyz[k] - g_z_y_0_yyz[k] * ab_y + g_z_y_0_yyyz[k];

                g_z_y_y_yzz[k] = g_z_0_0_yzz[k] - g_z_y_0_yzz[k] * ab_y + g_z_y_0_yyzz[k];

                g_z_y_y_zzz[k] = g_z_0_0_zzz[k] - g_z_y_0_zzz[k] * ab_y + g_z_y_0_yzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_z_y_z_xxx = cbuffer.data(pf_geom_11_off + 230 * ccomps * dcomps);

            auto g_z_y_z_xxy = cbuffer.data(pf_geom_11_off + 231 * ccomps * dcomps);

            auto g_z_y_z_xxz = cbuffer.data(pf_geom_11_off + 232 * ccomps * dcomps);

            auto g_z_y_z_xyy = cbuffer.data(pf_geom_11_off + 233 * ccomps * dcomps);

            auto g_z_y_z_xyz = cbuffer.data(pf_geom_11_off + 234 * ccomps * dcomps);

            auto g_z_y_z_xzz = cbuffer.data(pf_geom_11_off + 235 * ccomps * dcomps);

            auto g_z_y_z_yyy = cbuffer.data(pf_geom_11_off + 236 * ccomps * dcomps);

            auto g_z_y_z_yyz = cbuffer.data(pf_geom_11_off + 237 * ccomps * dcomps);

            auto g_z_y_z_yzz = cbuffer.data(pf_geom_11_off + 238 * ccomps * dcomps);

            auto g_z_y_z_zzz = cbuffer.data(pf_geom_11_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxx, g_0_y_0_xxy, g_0_y_0_xxz, g_0_y_0_xyy, g_0_y_0_xyz, g_0_y_0_xzz, g_0_y_0_yyy, g_0_y_0_yyz, g_0_y_0_yzz, g_0_y_0_zzz, g_z_y_0_xxx, g_z_y_0_xxxz, g_z_y_0_xxy, g_z_y_0_xxyz, g_z_y_0_xxz, g_z_y_0_xxzz, g_z_y_0_xyy, g_z_y_0_xyyz, g_z_y_0_xyz, g_z_y_0_xyzz, g_z_y_0_xzz, g_z_y_0_xzzz, g_z_y_0_yyy, g_z_y_0_yyyz, g_z_y_0_yyz, g_z_y_0_yyzz, g_z_y_0_yzz, g_z_y_0_yzzz, g_z_y_0_zzz, g_z_y_0_zzzz, g_z_y_z_xxx, g_z_y_z_xxy, g_z_y_z_xxz, g_z_y_z_xyy, g_z_y_z_xyz, g_z_y_z_xzz, g_z_y_z_yyy, g_z_y_z_yyz, g_z_y_z_yzz, g_z_y_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_z_xxx[k] = -g_0_y_0_xxx[k] - g_z_y_0_xxx[k] * ab_z + g_z_y_0_xxxz[k];

                g_z_y_z_xxy[k] = -g_0_y_0_xxy[k] - g_z_y_0_xxy[k] * ab_z + g_z_y_0_xxyz[k];

                g_z_y_z_xxz[k] = -g_0_y_0_xxz[k] - g_z_y_0_xxz[k] * ab_z + g_z_y_0_xxzz[k];

                g_z_y_z_xyy[k] = -g_0_y_0_xyy[k] - g_z_y_0_xyy[k] * ab_z + g_z_y_0_xyyz[k];

                g_z_y_z_xyz[k] = -g_0_y_0_xyz[k] - g_z_y_0_xyz[k] * ab_z + g_z_y_0_xyzz[k];

                g_z_y_z_xzz[k] = -g_0_y_0_xzz[k] - g_z_y_0_xzz[k] * ab_z + g_z_y_0_xzzz[k];

                g_z_y_z_yyy[k] = -g_0_y_0_yyy[k] - g_z_y_0_yyy[k] * ab_z + g_z_y_0_yyyz[k];

                g_z_y_z_yyz[k] = -g_0_y_0_yyz[k] - g_z_y_0_yyz[k] * ab_z + g_z_y_0_yyzz[k];

                g_z_y_z_yzz[k] = -g_0_y_0_yzz[k] - g_z_y_0_yzz[k] * ab_z + g_z_y_0_yzzz[k];

                g_z_y_z_zzz[k] = -g_0_y_0_zzz[k] - g_z_y_0_zzz[k] * ab_z + g_z_y_0_zzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_z_z_x_xxx = cbuffer.data(pf_geom_11_off + 240 * ccomps * dcomps);

            auto g_z_z_x_xxy = cbuffer.data(pf_geom_11_off + 241 * ccomps * dcomps);

            auto g_z_z_x_xxz = cbuffer.data(pf_geom_11_off + 242 * ccomps * dcomps);

            auto g_z_z_x_xyy = cbuffer.data(pf_geom_11_off + 243 * ccomps * dcomps);

            auto g_z_z_x_xyz = cbuffer.data(pf_geom_11_off + 244 * ccomps * dcomps);

            auto g_z_z_x_xzz = cbuffer.data(pf_geom_11_off + 245 * ccomps * dcomps);

            auto g_z_z_x_yyy = cbuffer.data(pf_geom_11_off + 246 * ccomps * dcomps);

            auto g_z_z_x_yyz = cbuffer.data(pf_geom_11_off + 247 * ccomps * dcomps);

            auto g_z_z_x_yzz = cbuffer.data(pf_geom_11_off + 248 * ccomps * dcomps);

            auto g_z_z_x_zzz = cbuffer.data(pf_geom_11_off + 249 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_0_xxx, g_z_z_0_xxxx, g_z_z_0_xxxy, g_z_z_0_xxxz, g_z_z_0_xxy, g_z_z_0_xxyy, g_z_z_0_xxyz, g_z_z_0_xxz, g_z_z_0_xxzz, g_z_z_0_xyy, g_z_z_0_xyyy, g_z_z_0_xyyz, g_z_z_0_xyz, g_z_z_0_xyzz, g_z_z_0_xzz, g_z_z_0_xzzz, g_z_z_0_yyy, g_z_z_0_yyz, g_z_z_0_yzz, g_z_z_0_zzz, g_z_z_x_xxx, g_z_z_x_xxy, g_z_z_x_xxz, g_z_z_x_xyy, g_z_z_x_xyz, g_z_z_x_xzz, g_z_z_x_yyy, g_z_z_x_yyz, g_z_z_x_yzz, g_z_z_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_x_xxx[k] = -g_z_z_0_xxx[k] * ab_x + g_z_z_0_xxxx[k];

                g_z_z_x_xxy[k] = -g_z_z_0_xxy[k] * ab_x + g_z_z_0_xxxy[k];

                g_z_z_x_xxz[k] = -g_z_z_0_xxz[k] * ab_x + g_z_z_0_xxxz[k];

                g_z_z_x_xyy[k] = -g_z_z_0_xyy[k] * ab_x + g_z_z_0_xxyy[k];

                g_z_z_x_xyz[k] = -g_z_z_0_xyz[k] * ab_x + g_z_z_0_xxyz[k];

                g_z_z_x_xzz[k] = -g_z_z_0_xzz[k] * ab_x + g_z_z_0_xxzz[k];

                g_z_z_x_yyy[k] = -g_z_z_0_yyy[k] * ab_x + g_z_z_0_xyyy[k];

                g_z_z_x_yyz[k] = -g_z_z_0_yyz[k] * ab_x + g_z_z_0_xyyz[k];

                g_z_z_x_yzz[k] = -g_z_z_0_yzz[k] * ab_x + g_z_z_0_xyzz[k];

                g_z_z_x_zzz[k] = -g_z_z_0_zzz[k] * ab_x + g_z_z_0_xzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_z_z_y_xxx = cbuffer.data(pf_geom_11_off + 250 * ccomps * dcomps);

            auto g_z_z_y_xxy = cbuffer.data(pf_geom_11_off + 251 * ccomps * dcomps);

            auto g_z_z_y_xxz = cbuffer.data(pf_geom_11_off + 252 * ccomps * dcomps);

            auto g_z_z_y_xyy = cbuffer.data(pf_geom_11_off + 253 * ccomps * dcomps);

            auto g_z_z_y_xyz = cbuffer.data(pf_geom_11_off + 254 * ccomps * dcomps);

            auto g_z_z_y_xzz = cbuffer.data(pf_geom_11_off + 255 * ccomps * dcomps);

            auto g_z_z_y_yyy = cbuffer.data(pf_geom_11_off + 256 * ccomps * dcomps);

            auto g_z_z_y_yyz = cbuffer.data(pf_geom_11_off + 257 * ccomps * dcomps);

            auto g_z_z_y_yzz = cbuffer.data(pf_geom_11_off + 258 * ccomps * dcomps);

            auto g_z_z_y_zzz = cbuffer.data(pf_geom_11_off + 259 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_0_xxx, g_z_z_0_xxxy, g_z_z_0_xxy, g_z_z_0_xxyy, g_z_z_0_xxyz, g_z_z_0_xxz, g_z_z_0_xyy, g_z_z_0_xyyy, g_z_z_0_xyyz, g_z_z_0_xyz, g_z_z_0_xyzz, g_z_z_0_xzz, g_z_z_0_yyy, g_z_z_0_yyyy, g_z_z_0_yyyz, g_z_z_0_yyz, g_z_z_0_yyzz, g_z_z_0_yzz, g_z_z_0_yzzz, g_z_z_0_zzz, g_z_z_y_xxx, g_z_z_y_xxy, g_z_z_y_xxz, g_z_z_y_xyy, g_z_z_y_xyz, g_z_z_y_xzz, g_z_z_y_yyy, g_z_z_y_yyz, g_z_z_y_yzz, g_z_z_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_y_xxx[k] = -g_z_z_0_xxx[k] * ab_y + g_z_z_0_xxxy[k];

                g_z_z_y_xxy[k] = -g_z_z_0_xxy[k] * ab_y + g_z_z_0_xxyy[k];

                g_z_z_y_xxz[k] = -g_z_z_0_xxz[k] * ab_y + g_z_z_0_xxyz[k];

                g_z_z_y_xyy[k] = -g_z_z_0_xyy[k] * ab_y + g_z_z_0_xyyy[k];

                g_z_z_y_xyz[k] = -g_z_z_0_xyz[k] * ab_y + g_z_z_0_xyyz[k];

                g_z_z_y_xzz[k] = -g_z_z_0_xzz[k] * ab_y + g_z_z_0_xyzz[k];

                g_z_z_y_yyy[k] = -g_z_z_0_yyy[k] * ab_y + g_z_z_0_yyyy[k];

                g_z_z_y_yyz[k] = -g_z_z_0_yyz[k] * ab_y + g_z_z_0_yyyz[k];

                g_z_z_y_yzz[k] = -g_z_z_0_yzz[k] * ab_y + g_z_z_0_yyzz[k];

                g_z_z_y_zzz[k] = -g_z_z_0_zzz[k] * ab_y + g_z_z_0_yzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_z_z_z_xxx = cbuffer.data(pf_geom_11_off + 260 * ccomps * dcomps);

            auto g_z_z_z_xxy = cbuffer.data(pf_geom_11_off + 261 * ccomps * dcomps);

            auto g_z_z_z_xxz = cbuffer.data(pf_geom_11_off + 262 * ccomps * dcomps);

            auto g_z_z_z_xyy = cbuffer.data(pf_geom_11_off + 263 * ccomps * dcomps);

            auto g_z_z_z_xyz = cbuffer.data(pf_geom_11_off + 264 * ccomps * dcomps);

            auto g_z_z_z_xzz = cbuffer.data(pf_geom_11_off + 265 * ccomps * dcomps);

            auto g_z_z_z_yyy = cbuffer.data(pf_geom_11_off + 266 * ccomps * dcomps);

            auto g_z_z_z_yyz = cbuffer.data(pf_geom_11_off + 267 * ccomps * dcomps);

            auto g_z_z_z_yzz = cbuffer.data(pf_geom_11_off + 268 * ccomps * dcomps);

            auto g_z_z_z_zzz = cbuffer.data(pf_geom_11_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxx, g_0_z_0_xxy, g_0_z_0_xxz, g_0_z_0_xyy, g_0_z_0_xyz, g_0_z_0_xzz, g_0_z_0_yyy, g_0_z_0_yyz, g_0_z_0_yzz, g_0_z_0_zzz, g_z_0_0_xxx, g_z_0_0_xxy, g_z_0_0_xxz, g_z_0_0_xyy, g_z_0_0_xyz, g_z_0_0_xzz, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yzz, g_z_0_0_zzz, g_z_z_0_xxx, g_z_z_0_xxxz, g_z_z_0_xxy, g_z_z_0_xxyz, g_z_z_0_xxz, g_z_z_0_xxzz, g_z_z_0_xyy, g_z_z_0_xyyz, g_z_z_0_xyz, g_z_z_0_xyzz, g_z_z_0_xzz, g_z_z_0_xzzz, g_z_z_0_yyy, g_z_z_0_yyyz, g_z_z_0_yyz, g_z_z_0_yyzz, g_z_z_0_yzz, g_z_z_0_yzzz, g_z_z_0_zzz, g_z_z_0_zzzz, g_z_z_z_xxx, g_z_z_z_xxy, g_z_z_z_xxz, g_z_z_z_xyy, g_z_z_z_xyz, g_z_z_z_xzz, g_z_z_z_yyy, g_z_z_z_yyz, g_z_z_z_yzz, g_z_z_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_z_xxx[k] = -g_0_z_0_xxx[k] + g_z_0_0_xxx[k] - g_z_z_0_xxx[k] * ab_z + g_z_z_0_xxxz[k];

                g_z_z_z_xxy[k] = -g_0_z_0_xxy[k] + g_z_0_0_xxy[k] - g_z_z_0_xxy[k] * ab_z + g_z_z_0_xxyz[k];

                g_z_z_z_xxz[k] = -g_0_z_0_xxz[k] + g_z_0_0_xxz[k] - g_z_z_0_xxz[k] * ab_z + g_z_z_0_xxzz[k];

                g_z_z_z_xyy[k] = -g_0_z_0_xyy[k] + g_z_0_0_xyy[k] - g_z_z_0_xyy[k] * ab_z + g_z_z_0_xyyz[k];

                g_z_z_z_xyz[k] = -g_0_z_0_xyz[k] + g_z_0_0_xyz[k] - g_z_z_0_xyz[k] * ab_z + g_z_z_0_xyzz[k];

                g_z_z_z_xzz[k] = -g_0_z_0_xzz[k] + g_z_0_0_xzz[k] - g_z_z_0_xzz[k] * ab_z + g_z_z_0_xzzz[k];

                g_z_z_z_yyy[k] = -g_0_z_0_yyy[k] + g_z_0_0_yyy[k] - g_z_z_0_yyy[k] * ab_z + g_z_z_0_yyyz[k];

                g_z_z_z_yyz[k] = -g_0_z_0_yyz[k] + g_z_0_0_yyz[k] - g_z_z_0_yyz[k] * ab_z + g_z_z_0_yyzz[k];

                g_z_z_z_yzz[k] = -g_0_z_0_yzz[k] + g_z_0_0_yzz[k] - g_z_z_0_yzz[k] * ab_z + g_z_z_0_yzzz[k];

                g_z_z_z_zzz[k] = -g_0_z_0_zzz[k] + g_z_0_0_zzz[k] - g_z_z_0_zzz[k] * ab_z + g_z_z_0_zzzz[k];
            }
        }
    }
}

} // erirec namespace

