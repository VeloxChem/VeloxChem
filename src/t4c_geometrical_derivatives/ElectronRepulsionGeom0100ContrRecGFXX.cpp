#include "ElectronRepulsionGeom0100ContrRecGFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_gfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_gfxx,
                                            const size_t idx_ffxx,
                                            const size_t idx_geom_01_ffxx,
                                            const size_t idx_geom_01_fgxx,
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
            /// Set up components of auxilary buffer : FFSS

            const auto ff_off = idx_ffxx + i * dcomps + j;

            auto g_xxx_xxx = cbuffer.data(ff_off + 0 * ccomps * dcomps);

            auto g_xxx_xxy = cbuffer.data(ff_off + 1 * ccomps * dcomps);

            auto g_xxx_xxz = cbuffer.data(ff_off + 2 * ccomps * dcomps);

            auto g_xxx_xyy = cbuffer.data(ff_off + 3 * ccomps * dcomps);

            auto g_xxx_xyz = cbuffer.data(ff_off + 4 * ccomps * dcomps);

            auto g_xxx_xzz = cbuffer.data(ff_off + 5 * ccomps * dcomps);

            auto g_xxx_yyy = cbuffer.data(ff_off + 6 * ccomps * dcomps);

            auto g_xxx_yyz = cbuffer.data(ff_off + 7 * ccomps * dcomps);

            auto g_xxx_yzz = cbuffer.data(ff_off + 8 * ccomps * dcomps);

            auto g_xxx_zzz = cbuffer.data(ff_off + 9 * ccomps * dcomps);

            auto g_xxy_xxx = cbuffer.data(ff_off + 10 * ccomps * dcomps);

            auto g_xxy_xxy = cbuffer.data(ff_off + 11 * ccomps * dcomps);

            auto g_xxy_xxz = cbuffer.data(ff_off + 12 * ccomps * dcomps);

            auto g_xxy_xyy = cbuffer.data(ff_off + 13 * ccomps * dcomps);

            auto g_xxy_xyz = cbuffer.data(ff_off + 14 * ccomps * dcomps);

            auto g_xxy_xzz = cbuffer.data(ff_off + 15 * ccomps * dcomps);

            auto g_xxy_yyy = cbuffer.data(ff_off + 16 * ccomps * dcomps);

            auto g_xxy_yyz = cbuffer.data(ff_off + 17 * ccomps * dcomps);

            auto g_xxy_yzz = cbuffer.data(ff_off + 18 * ccomps * dcomps);

            auto g_xxy_zzz = cbuffer.data(ff_off + 19 * ccomps * dcomps);

            auto g_xxz_xxx = cbuffer.data(ff_off + 20 * ccomps * dcomps);

            auto g_xxz_xxy = cbuffer.data(ff_off + 21 * ccomps * dcomps);

            auto g_xxz_xxz = cbuffer.data(ff_off + 22 * ccomps * dcomps);

            auto g_xxz_xyy = cbuffer.data(ff_off + 23 * ccomps * dcomps);

            auto g_xxz_xyz = cbuffer.data(ff_off + 24 * ccomps * dcomps);

            auto g_xxz_xzz = cbuffer.data(ff_off + 25 * ccomps * dcomps);

            auto g_xxz_yyy = cbuffer.data(ff_off + 26 * ccomps * dcomps);

            auto g_xxz_yyz = cbuffer.data(ff_off + 27 * ccomps * dcomps);

            auto g_xxz_yzz = cbuffer.data(ff_off + 28 * ccomps * dcomps);

            auto g_xxz_zzz = cbuffer.data(ff_off + 29 * ccomps * dcomps);

            auto g_xyy_xxx = cbuffer.data(ff_off + 30 * ccomps * dcomps);

            auto g_xyy_xxy = cbuffer.data(ff_off + 31 * ccomps * dcomps);

            auto g_xyy_xxz = cbuffer.data(ff_off + 32 * ccomps * dcomps);

            auto g_xyy_xyy = cbuffer.data(ff_off + 33 * ccomps * dcomps);

            auto g_xyy_xyz = cbuffer.data(ff_off + 34 * ccomps * dcomps);

            auto g_xyy_xzz = cbuffer.data(ff_off + 35 * ccomps * dcomps);

            auto g_xyy_yyy = cbuffer.data(ff_off + 36 * ccomps * dcomps);

            auto g_xyy_yyz = cbuffer.data(ff_off + 37 * ccomps * dcomps);

            auto g_xyy_yzz = cbuffer.data(ff_off + 38 * ccomps * dcomps);

            auto g_xyy_zzz = cbuffer.data(ff_off + 39 * ccomps * dcomps);

            auto g_xyz_xxx = cbuffer.data(ff_off + 40 * ccomps * dcomps);

            auto g_xyz_xxy = cbuffer.data(ff_off + 41 * ccomps * dcomps);

            auto g_xyz_xxz = cbuffer.data(ff_off + 42 * ccomps * dcomps);

            auto g_xyz_xyy = cbuffer.data(ff_off + 43 * ccomps * dcomps);

            auto g_xyz_xyz = cbuffer.data(ff_off + 44 * ccomps * dcomps);

            auto g_xyz_xzz = cbuffer.data(ff_off + 45 * ccomps * dcomps);

            auto g_xyz_yyy = cbuffer.data(ff_off + 46 * ccomps * dcomps);

            auto g_xyz_yyz = cbuffer.data(ff_off + 47 * ccomps * dcomps);

            auto g_xyz_yzz = cbuffer.data(ff_off + 48 * ccomps * dcomps);

            auto g_xyz_zzz = cbuffer.data(ff_off + 49 * ccomps * dcomps);

            auto g_xzz_xxx = cbuffer.data(ff_off + 50 * ccomps * dcomps);

            auto g_xzz_xxy = cbuffer.data(ff_off + 51 * ccomps * dcomps);

            auto g_xzz_xxz = cbuffer.data(ff_off + 52 * ccomps * dcomps);

            auto g_xzz_xyy = cbuffer.data(ff_off + 53 * ccomps * dcomps);

            auto g_xzz_xyz = cbuffer.data(ff_off + 54 * ccomps * dcomps);

            auto g_xzz_xzz = cbuffer.data(ff_off + 55 * ccomps * dcomps);

            auto g_xzz_yyy = cbuffer.data(ff_off + 56 * ccomps * dcomps);

            auto g_xzz_yyz = cbuffer.data(ff_off + 57 * ccomps * dcomps);

            auto g_xzz_yzz = cbuffer.data(ff_off + 58 * ccomps * dcomps);

            auto g_xzz_zzz = cbuffer.data(ff_off + 59 * ccomps * dcomps);

            auto g_yyy_xxx = cbuffer.data(ff_off + 60 * ccomps * dcomps);

            auto g_yyy_xxy = cbuffer.data(ff_off + 61 * ccomps * dcomps);

            auto g_yyy_xxz = cbuffer.data(ff_off + 62 * ccomps * dcomps);

            auto g_yyy_xyy = cbuffer.data(ff_off + 63 * ccomps * dcomps);

            auto g_yyy_xyz = cbuffer.data(ff_off + 64 * ccomps * dcomps);

            auto g_yyy_xzz = cbuffer.data(ff_off + 65 * ccomps * dcomps);

            auto g_yyy_yyy = cbuffer.data(ff_off + 66 * ccomps * dcomps);

            auto g_yyy_yyz = cbuffer.data(ff_off + 67 * ccomps * dcomps);

            auto g_yyy_yzz = cbuffer.data(ff_off + 68 * ccomps * dcomps);

            auto g_yyy_zzz = cbuffer.data(ff_off + 69 * ccomps * dcomps);

            auto g_yyz_xxx = cbuffer.data(ff_off + 70 * ccomps * dcomps);

            auto g_yyz_xxy = cbuffer.data(ff_off + 71 * ccomps * dcomps);

            auto g_yyz_xxz = cbuffer.data(ff_off + 72 * ccomps * dcomps);

            auto g_yyz_xyy = cbuffer.data(ff_off + 73 * ccomps * dcomps);

            auto g_yyz_xyz = cbuffer.data(ff_off + 74 * ccomps * dcomps);

            auto g_yyz_xzz = cbuffer.data(ff_off + 75 * ccomps * dcomps);

            auto g_yyz_yyy = cbuffer.data(ff_off + 76 * ccomps * dcomps);

            auto g_yyz_yyz = cbuffer.data(ff_off + 77 * ccomps * dcomps);

            auto g_yyz_yzz = cbuffer.data(ff_off + 78 * ccomps * dcomps);

            auto g_yyz_zzz = cbuffer.data(ff_off + 79 * ccomps * dcomps);

            auto g_yzz_xxx = cbuffer.data(ff_off + 80 * ccomps * dcomps);

            auto g_yzz_xxy = cbuffer.data(ff_off + 81 * ccomps * dcomps);

            auto g_yzz_xxz = cbuffer.data(ff_off + 82 * ccomps * dcomps);

            auto g_yzz_xyy = cbuffer.data(ff_off + 83 * ccomps * dcomps);

            auto g_yzz_xyz = cbuffer.data(ff_off + 84 * ccomps * dcomps);

            auto g_yzz_xzz = cbuffer.data(ff_off + 85 * ccomps * dcomps);

            auto g_yzz_yyy = cbuffer.data(ff_off + 86 * ccomps * dcomps);

            auto g_yzz_yyz = cbuffer.data(ff_off + 87 * ccomps * dcomps);

            auto g_yzz_yzz = cbuffer.data(ff_off + 88 * ccomps * dcomps);

            auto g_yzz_zzz = cbuffer.data(ff_off + 89 * ccomps * dcomps);

            auto g_zzz_xxx = cbuffer.data(ff_off + 90 * ccomps * dcomps);

            auto g_zzz_xxy = cbuffer.data(ff_off + 91 * ccomps * dcomps);

            auto g_zzz_xxz = cbuffer.data(ff_off + 92 * ccomps * dcomps);

            auto g_zzz_xyy = cbuffer.data(ff_off + 93 * ccomps * dcomps);

            auto g_zzz_xyz = cbuffer.data(ff_off + 94 * ccomps * dcomps);

            auto g_zzz_xzz = cbuffer.data(ff_off + 95 * ccomps * dcomps);

            auto g_zzz_yyy = cbuffer.data(ff_off + 96 * ccomps * dcomps);

            auto g_zzz_yyz = cbuffer.data(ff_off + 97 * ccomps * dcomps);

            auto g_zzz_yzz = cbuffer.data(ff_off + 98 * ccomps * dcomps);

            auto g_zzz_zzz = cbuffer.data(ff_off + 99 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FFSS

            const auto ff_geom_01_off = idx_geom_01_ffxx + i * dcomps + j;

            auto g_0_x_xxx_xxx = cbuffer.data(ff_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xxy = cbuffer.data(ff_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xxz = cbuffer.data(ff_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_xyy = cbuffer.data(ff_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_xyz = cbuffer.data(ff_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_xzz = cbuffer.data(ff_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxx_yyy = cbuffer.data(ff_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxx_yyz = cbuffer.data(ff_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxx_yzz = cbuffer.data(ff_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxx_zzz = cbuffer.data(ff_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxy_xxx = cbuffer.data(ff_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxy_xxy = cbuffer.data(ff_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxy_xxz = cbuffer.data(ff_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxy_xyy = cbuffer.data(ff_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxy_xyz = cbuffer.data(ff_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxy_xzz = cbuffer.data(ff_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxy_yyy = cbuffer.data(ff_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxy_yyz = cbuffer.data(ff_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxy_yzz = cbuffer.data(ff_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxy_zzz = cbuffer.data(ff_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxz_xxx = cbuffer.data(ff_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxz_xxy = cbuffer.data(ff_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxz_xxz = cbuffer.data(ff_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxz_xyy = cbuffer.data(ff_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxz_xyz = cbuffer.data(ff_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxz_xzz = cbuffer.data(ff_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxz_yyy = cbuffer.data(ff_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxz_yyz = cbuffer.data(ff_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxz_yzz = cbuffer.data(ff_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxz_zzz = cbuffer.data(ff_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xyy_xxx = cbuffer.data(ff_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xyy_xxy = cbuffer.data(ff_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xyy_xxz = cbuffer.data(ff_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xyy_xyy = cbuffer.data(ff_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xyy_xyz = cbuffer.data(ff_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xyy_xzz = cbuffer.data(ff_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xyy_yyy = cbuffer.data(ff_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xyy_yyz = cbuffer.data(ff_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xyy_yzz = cbuffer.data(ff_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xyy_zzz = cbuffer.data(ff_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xyz_xxx = cbuffer.data(ff_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xyz_xxy = cbuffer.data(ff_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xyz_xxz = cbuffer.data(ff_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xyz_xyy = cbuffer.data(ff_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xyz_xyz = cbuffer.data(ff_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xyz_xzz = cbuffer.data(ff_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xyz_yyy = cbuffer.data(ff_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xyz_yyz = cbuffer.data(ff_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xyz_yzz = cbuffer.data(ff_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xyz_zzz = cbuffer.data(ff_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xzz_xxx = cbuffer.data(ff_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xzz_xxy = cbuffer.data(ff_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xzz_xxz = cbuffer.data(ff_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xzz_xyy = cbuffer.data(ff_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xzz_xyz = cbuffer.data(ff_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xzz_xzz = cbuffer.data(ff_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xzz_yyy = cbuffer.data(ff_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xzz_yyz = cbuffer.data(ff_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xzz_yzz = cbuffer.data(ff_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xzz_zzz = cbuffer.data(ff_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_yyy_xxx = cbuffer.data(ff_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_yyy_xxy = cbuffer.data(ff_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_yyy_xxz = cbuffer.data(ff_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_yyy_xyy = cbuffer.data(ff_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yyy_xyz = cbuffer.data(ff_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yyy_xzz = cbuffer.data(ff_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_yyy_yyy = cbuffer.data(ff_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yyy_yyz = cbuffer.data(ff_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yyy_yzz = cbuffer.data(ff_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yyy_zzz = cbuffer.data(ff_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yyz_xxx = cbuffer.data(ff_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yyz_xxy = cbuffer.data(ff_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_yyz_xxz = cbuffer.data(ff_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yyz_xyy = cbuffer.data(ff_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yyz_xyz = cbuffer.data(ff_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_yyz_xzz = cbuffer.data(ff_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_yyz_yyy = cbuffer.data(ff_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_yyz_yyz = cbuffer.data(ff_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_yyz_yzz = cbuffer.data(ff_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_yyz_zzz = cbuffer.data(ff_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_yzz_xxx = cbuffer.data(ff_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_yzz_xxy = cbuffer.data(ff_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_yzz_xxz = cbuffer.data(ff_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_yzz_xyy = cbuffer.data(ff_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_yzz_xyz = cbuffer.data(ff_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_yzz_xzz = cbuffer.data(ff_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_yzz_yyy = cbuffer.data(ff_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_yzz_yyz = cbuffer.data(ff_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_yzz_yzz = cbuffer.data(ff_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_yzz_zzz = cbuffer.data(ff_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_zzz_xxx = cbuffer.data(ff_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_zzz_xxy = cbuffer.data(ff_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_zzz_xxz = cbuffer.data(ff_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_zzz_xyy = cbuffer.data(ff_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_zzz_xyz = cbuffer.data(ff_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_zzz_xzz = cbuffer.data(ff_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_zzz_yyy = cbuffer.data(ff_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_zzz_yyz = cbuffer.data(ff_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_zzz_yzz = cbuffer.data(ff_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_zzz_zzz = cbuffer.data(ff_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xxx_xxx = cbuffer.data(ff_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xxx_xxy = cbuffer.data(ff_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_xxx_xxz = cbuffer.data(ff_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xxx_xyy = cbuffer.data(ff_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xxx_xyz = cbuffer.data(ff_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_xxx_xzz = cbuffer.data(ff_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xxx_yyy = cbuffer.data(ff_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xxx_yyz = cbuffer.data(ff_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_xxx_yzz = cbuffer.data(ff_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xxx_zzz = cbuffer.data(ff_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xxy_xxx = cbuffer.data(ff_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_xxy_xxy = cbuffer.data(ff_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xxy_xxz = cbuffer.data(ff_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xxy_xyy = cbuffer.data(ff_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_xxy_xyz = cbuffer.data(ff_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xxy_xzz = cbuffer.data(ff_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xxy_yyy = cbuffer.data(ff_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_xxy_yyz = cbuffer.data(ff_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xxy_yzz = cbuffer.data(ff_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xxy_zzz = cbuffer.data(ff_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_xxz_xxx = cbuffer.data(ff_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xxz_xxy = cbuffer.data(ff_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xxz_xxz = cbuffer.data(ff_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_xxz_xyy = cbuffer.data(ff_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xxz_xyz = cbuffer.data(ff_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xxz_xzz = cbuffer.data(ff_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xxz_yyy = cbuffer.data(ff_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xxz_yyz = cbuffer.data(ff_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xxz_yzz = cbuffer.data(ff_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xxz_zzz = cbuffer.data(ff_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xyy_xxx = cbuffer.data(ff_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xyy_xxy = cbuffer.data(ff_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xyy_xxz = cbuffer.data(ff_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xyy_xyy = cbuffer.data(ff_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xyy_xyz = cbuffer.data(ff_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xyy_xzz = cbuffer.data(ff_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xyy_yyy = cbuffer.data(ff_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xyy_yyz = cbuffer.data(ff_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_xyy_yzz = cbuffer.data(ff_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xyy_zzz = cbuffer.data(ff_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xyz_xxx = cbuffer.data(ff_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xyz_xxy = cbuffer.data(ff_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xyz_xxz = cbuffer.data(ff_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xyz_xyy = cbuffer.data(ff_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_xyz_xyz = cbuffer.data(ff_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xyz_xzz = cbuffer.data(ff_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xyz_yyy = cbuffer.data(ff_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_xyz_yyz = cbuffer.data(ff_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xyz_yzz = cbuffer.data(ff_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xyz_zzz = cbuffer.data(ff_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_xzz_xxx = cbuffer.data(ff_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xzz_xxy = cbuffer.data(ff_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xzz_xxz = cbuffer.data(ff_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xzz_xyy = cbuffer.data(ff_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xzz_xyz = cbuffer.data(ff_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xzz_xzz = cbuffer.data(ff_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xzz_yyy = cbuffer.data(ff_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xzz_yyz = cbuffer.data(ff_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xzz_yzz = cbuffer.data(ff_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xzz_zzz = cbuffer.data(ff_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_yyy_xxx = cbuffer.data(ff_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_yyy_xxy = cbuffer.data(ff_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_yyy_xxz = cbuffer.data(ff_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_yyy_xyy = cbuffer.data(ff_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_yyy_xyz = cbuffer.data(ff_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_yyy_xzz = cbuffer.data(ff_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_yyy_yyy = cbuffer.data(ff_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_yyy_yyz = cbuffer.data(ff_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_yyy_yzz = cbuffer.data(ff_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_yyy_zzz = cbuffer.data(ff_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_yyz_xxx = cbuffer.data(ff_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_yyz_xxy = cbuffer.data(ff_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_yyz_xxz = cbuffer.data(ff_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_yyz_xyy = cbuffer.data(ff_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_yyz_xyz = cbuffer.data(ff_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_yyz_xzz = cbuffer.data(ff_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_yyz_yyy = cbuffer.data(ff_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_yyz_yyz = cbuffer.data(ff_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_yyz_yzz = cbuffer.data(ff_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_yyz_zzz = cbuffer.data(ff_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_yzz_xxx = cbuffer.data(ff_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_yzz_xxy = cbuffer.data(ff_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_yzz_xxz = cbuffer.data(ff_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_yzz_xyy = cbuffer.data(ff_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_yzz_xyz = cbuffer.data(ff_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_yzz_xzz = cbuffer.data(ff_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_yzz_yyy = cbuffer.data(ff_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_yzz_yyz = cbuffer.data(ff_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_yzz_yzz = cbuffer.data(ff_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_yzz_zzz = cbuffer.data(ff_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_zzz_xxx = cbuffer.data(ff_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_zzz_xxy = cbuffer.data(ff_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_zzz_xxz = cbuffer.data(ff_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_zzz_xyy = cbuffer.data(ff_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_zzz_xyz = cbuffer.data(ff_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_zzz_xzz = cbuffer.data(ff_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_zzz_yyy = cbuffer.data(ff_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_zzz_yyz = cbuffer.data(ff_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_zzz_yzz = cbuffer.data(ff_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_zzz_zzz = cbuffer.data(ff_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_xxx_xxx = cbuffer.data(ff_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_xxx_xxy = cbuffer.data(ff_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_xxx_xxz = cbuffer.data(ff_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_xxx_xyy = cbuffer.data(ff_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_z_xxx_xyz = cbuffer.data(ff_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_xxx_xzz = cbuffer.data(ff_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_xxx_yyy = cbuffer.data(ff_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_xxx_yyz = cbuffer.data(ff_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_xxx_yzz = cbuffer.data(ff_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_xxx_zzz = cbuffer.data(ff_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_z_xxy_xxx = cbuffer.data(ff_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_xxy_xxy = cbuffer.data(ff_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_xxy_xxz = cbuffer.data(ff_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_xxy_xyy = cbuffer.data(ff_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_xxy_xyz = cbuffer.data(ff_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_xxy_xzz = cbuffer.data(ff_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_xxy_yyy = cbuffer.data(ff_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xxy_yyz = cbuffer.data(ff_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xxy_yzz = cbuffer.data(ff_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_xxy_zzz = cbuffer.data(ff_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xxz_xxx = cbuffer.data(ff_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xxz_xxy = cbuffer.data(ff_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_xxz_xxz = cbuffer.data(ff_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xxz_xyy = cbuffer.data(ff_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xxz_xyz = cbuffer.data(ff_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_xxz_xzz = cbuffer.data(ff_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_xxz_yyy = cbuffer.data(ff_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_xxz_yyz = cbuffer.data(ff_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_xxz_yzz = cbuffer.data(ff_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_xxz_zzz = cbuffer.data(ff_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_xyy_xxx = cbuffer.data(ff_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_xyy_xxy = cbuffer.data(ff_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_xyy_xxz = cbuffer.data(ff_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_xyy_xyy = cbuffer.data(ff_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_xyy_xyz = cbuffer.data(ff_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_xyy_xzz = cbuffer.data(ff_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_xyy_yyy = cbuffer.data(ff_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_xyy_yyz = cbuffer.data(ff_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_xyy_yzz = cbuffer.data(ff_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_xyy_zzz = cbuffer.data(ff_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_xyz_xxx = cbuffer.data(ff_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_xyz_xxy = cbuffer.data(ff_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_xyz_xxz = cbuffer.data(ff_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_xyz_xyy = cbuffer.data(ff_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_xyz_xyz = cbuffer.data(ff_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_xyz_xzz = cbuffer.data(ff_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_xyz_yyy = cbuffer.data(ff_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_xyz_yyz = cbuffer.data(ff_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_xyz_yzz = cbuffer.data(ff_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_xyz_zzz = cbuffer.data(ff_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_xzz_xxx = cbuffer.data(ff_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_xzz_xxy = cbuffer.data(ff_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_xzz_xxz = cbuffer.data(ff_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_xzz_xyy = cbuffer.data(ff_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_xzz_xyz = cbuffer.data(ff_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_xzz_xzz = cbuffer.data(ff_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_xzz_yyy = cbuffer.data(ff_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_xzz_yyz = cbuffer.data(ff_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_xzz_yzz = cbuffer.data(ff_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_xzz_zzz = cbuffer.data(ff_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_yyy_xxx = cbuffer.data(ff_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_yyy_xxy = cbuffer.data(ff_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_yyy_xxz = cbuffer.data(ff_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_yyy_xyy = cbuffer.data(ff_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_yyy_xyz = cbuffer.data(ff_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_yyy_xzz = cbuffer.data(ff_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_yyy_yyy = cbuffer.data(ff_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_yyy_yyz = cbuffer.data(ff_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_yyy_yzz = cbuffer.data(ff_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_yyy_zzz = cbuffer.data(ff_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_z_yyz_xxx = cbuffer.data(ff_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_yyz_xxy = cbuffer.data(ff_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_yyz_xxz = cbuffer.data(ff_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_z_yyz_xyy = cbuffer.data(ff_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_yyz_xyz = cbuffer.data(ff_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_yyz_xzz = cbuffer.data(ff_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_z_yyz_yyy = cbuffer.data(ff_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_yyz_yyz = cbuffer.data(ff_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_yyz_yzz = cbuffer.data(ff_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_z_yyz_zzz = cbuffer.data(ff_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_yzz_xxx = cbuffer.data(ff_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_yzz_xxy = cbuffer.data(ff_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_z_yzz_xxz = cbuffer.data(ff_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_yzz_xyy = cbuffer.data(ff_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_yzz_xyz = cbuffer.data(ff_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_z_yzz_xzz = cbuffer.data(ff_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_yzz_yyy = cbuffer.data(ff_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_yzz_yyz = cbuffer.data(ff_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_z_yzz_yzz = cbuffer.data(ff_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_yzz_zzz = cbuffer.data(ff_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_zzz_xxx = cbuffer.data(ff_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_z_zzz_xxy = cbuffer.data(ff_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_zzz_xxz = cbuffer.data(ff_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_zzz_xyy = cbuffer.data(ff_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_z_zzz_xyz = cbuffer.data(ff_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_zzz_xzz = cbuffer.data(ff_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_zzz_yyy = cbuffer.data(ff_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_z_zzz_yyz = cbuffer.data(ff_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_zzz_yzz = cbuffer.data(ff_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_zzz_zzz = cbuffer.data(ff_geom_01_off + 299 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FGSS

            const auto fg_geom_01_off = idx_geom_01_fgxx + i * dcomps + j;

            auto g_0_x_xxx_xxxx = cbuffer.data(fg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xxxy = cbuffer.data(fg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xxxz = cbuffer.data(fg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_xxyy = cbuffer.data(fg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_xxyz = cbuffer.data(fg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_xxzz = cbuffer.data(fg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxx_xyyy = cbuffer.data(fg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxx_xyyz = cbuffer.data(fg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxx_xyzz = cbuffer.data(fg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxx_xzzz = cbuffer.data(fg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxx_yyyy = cbuffer.data(fg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxx_yyyz = cbuffer.data(fg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxx_yyzz = cbuffer.data(fg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxx_yzzz = cbuffer.data(fg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxx_zzzz = cbuffer.data(fg_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxy_xxxx = cbuffer.data(fg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxy_xxxy = cbuffer.data(fg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxy_xxxz = cbuffer.data(fg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxy_xxyy = cbuffer.data(fg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxy_xxyz = cbuffer.data(fg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxy_xxzz = cbuffer.data(fg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxy_xyyy = cbuffer.data(fg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxy_xyyz = cbuffer.data(fg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxy_xyzz = cbuffer.data(fg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxy_xzzz = cbuffer.data(fg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxy_yyyy = cbuffer.data(fg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxy_yyyz = cbuffer.data(fg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxy_yyzz = cbuffer.data(fg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxy_yzzz = cbuffer.data(fg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxy_zzzz = cbuffer.data(fg_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxz_xxxx = cbuffer.data(fg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxz_xxxy = cbuffer.data(fg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxz_xxxz = cbuffer.data(fg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxz_xxyy = cbuffer.data(fg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxz_xxyz = cbuffer.data(fg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxz_xxzz = cbuffer.data(fg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxz_xyyy = cbuffer.data(fg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxz_xyyz = cbuffer.data(fg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxz_xyzz = cbuffer.data(fg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxz_xzzz = cbuffer.data(fg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxz_yyyy = cbuffer.data(fg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxz_yyyz = cbuffer.data(fg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxz_yyzz = cbuffer.data(fg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxz_yzzz = cbuffer.data(fg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxz_zzzz = cbuffer.data(fg_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xyy_xxxx = cbuffer.data(fg_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xyy_xxxy = cbuffer.data(fg_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xyy_xxxz = cbuffer.data(fg_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xyy_xxyy = cbuffer.data(fg_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xyy_xxyz = cbuffer.data(fg_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xyy_xxzz = cbuffer.data(fg_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xyy_xyyy = cbuffer.data(fg_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xyy_xyyz = cbuffer.data(fg_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xyy_xyzz = cbuffer.data(fg_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xyy_xzzz = cbuffer.data(fg_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xyy_yyyy = cbuffer.data(fg_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xyy_yyyz = cbuffer.data(fg_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xyy_yyzz = cbuffer.data(fg_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xyy_yzzz = cbuffer.data(fg_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xyy_zzzz = cbuffer.data(fg_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xyz_xxxx = cbuffer.data(fg_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xyz_xxxy = cbuffer.data(fg_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xyz_xxxz = cbuffer.data(fg_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xyz_xxyy = cbuffer.data(fg_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyz_xxyz = cbuffer.data(fg_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyz_xxzz = cbuffer.data(fg_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xyz_xyyy = cbuffer.data(fg_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyz_xyyz = cbuffer.data(fg_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyz_xyzz = cbuffer.data(fg_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xyz_xzzz = cbuffer.data(fg_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xyz_yyyy = cbuffer.data(fg_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyz_yyyz = cbuffer.data(fg_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xyz_yyzz = cbuffer.data(fg_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyz_yzzz = cbuffer.data(fg_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyz_zzzz = cbuffer.data(fg_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xzz_xxxx = cbuffer.data(fg_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xzz_xxxy = cbuffer.data(fg_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xzz_xxxz = cbuffer.data(fg_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xzz_xxyy = cbuffer.data(fg_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xzz_xxyz = cbuffer.data(fg_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xzz_xxzz = cbuffer.data(fg_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xzz_xyyy = cbuffer.data(fg_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xzz_xyyz = cbuffer.data(fg_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xzz_xyzz = cbuffer.data(fg_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xzz_xzzz = cbuffer.data(fg_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xzz_yyyy = cbuffer.data(fg_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xzz_yyyz = cbuffer.data(fg_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xzz_yyzz = cbuffer.data(fg_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xzz_yzzz = cbuffer.data(fg_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xzz_zzzz = cbuffer.data(fg_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_yyy_xxxx = cbuffer.data(fg_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_yyy_xxxy = cbuffer.data(fg_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_yyy_xxxz = cbuffer.data(fg_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_yyy_xxyy = cbuffer.data(fg_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_yyy_xxyz = cbuffer.data(fg_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_yyy_xxzz = cbuffer.data(fg_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_yyy_xyyy = cbuffer.data(fg_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_yyy_xyyz = cbuffer.data(fg_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_yyy_xyzz = cbuffer.data(fg_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_yyy_xzzz = cbuffer.data(fg_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_yyy_yyyy = cbuffer.data(fg_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yyy_yyyz = cbuffer.data(fg_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_yyy_yyzz = cbuffer.data(fg_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yyy_yzzz = cbuffer.data(fg_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yyy_zzzz = cbuffer.data(fg_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_yyz_xxxx = cbuffer.data(fg_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_yyz_xxxy = cbuffer.data(fg_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_yyz_xxxz = cbuffer.data(fg_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_yyz_xxyy = cbuffer.data(fg_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_yyz_xxyz = cbuffer.data(fg_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_yyz_xxzz = cbuffer.data(fg_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_yyz_xyyy = cbuffer.data(fg_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_yyz_xyyz = cbuffer.data(fg_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_yyz_xyzz = cbuffer.data(fg_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_yyz_xzzz = cbuffer.data(fg_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_yyz_yyyy = cbuffer.data(fg_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_yyz_yyyz = cbuffer.data(fg_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_yyz_yyzz = cbuffer.data(fg_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_yyz_yzzz = cbuffer.data(fg_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_yyz_zzzz = cbuffer.data(fg_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_yzz_xxxx = cbuffer.data(fg_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_yzz_xxxy = cbuffer.data(fg_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_yzz_xxxz = cbuffer.data(fg_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_yzz_xxyy = cbuffer.data(fg_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_yzz_xxyz = cbuffer.data(fg_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_yzz_xxzz = cbuffer.data(fg_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yzz_xyyy = cbuffer.data(fg_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yzz_xyyz = cbuffer.data(fg_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yzz_xyzz = cbuffer.data(fg_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yzz_xzzz = cbuffer.data(fg_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yzz_yyyy = cbuffer.data(fg_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yzz_yyyz = cbuffer.data(fg_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yzz_yyzz = cbuffer.data(fg_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yzz_yzzz = cbuffer.data(fg_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yzz_zzzz = cbuffer.data(fg_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_zzz_xxxx = cbuffer.data(fg_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_zzz_xxxy = cbuffer.data(fg_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_zzz_xxxz = cbuffer.data(fg_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_zzz_xxyy = cbuffer.data(fg_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_zzz_xxyz = cbuffer.data(fg_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_zzz_xxzz = cbuffer.data(fg_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_zzz_xyyy = cbuffer.data(fg_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_zzz_xyyz = cbuffer.data(fg_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_zzz_xyzz = cbuffer.data(fg_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_zzz_xzzz = cbuffer.data(fg_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_zzz_yyyy = cbuffer.data(fg_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_zzz_yyyz = cbuffer.data(fg_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_zzz_yyzz = cbuffer.data(fg_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_zzz_yzzz = cbuffer.data(fg_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_zzz_zzzz = cbuffer.data(fg_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_xxx_xxxx = cbuffer.data(fg_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xxx_xxxy = cbuffer.data(fg_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xxx_xxxz = cbuffer.data(fg_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xxx_xxyy = cbuffer.data(fg_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xxx_xxyz = cbuffer.data(fg_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xxx_xxzz = cbuffer.data(fg_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xxx_xyyy = cbuffer.data(fg_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xxx_xyyz = cbuffer.data(fg_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xxx_xyzz = cbuffer.data(fg_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xxx_xzzz = cbuffer.data(fg_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_xxx_yyyy = cbuffer.data(fg_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xxx_yyyz = cbuffer.data(fg_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_xxx_yyzz = cbuffer.data(fg_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xxx_yzzz = cbuffer.data(fg_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xxx_zzzz = cbuffer.data(fg_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_xxy_xxxx = cbuffer.data(fg_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xxy_xxxy = cbuffer.data(fg_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xxy_xxxz = cbuffer.data(fg_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_xxy_xxyy = cbuffer.data(fg_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxy_xxyz = cbuffer.data(fg_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxy_xxzz = cbuffer.data(fg_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxy_xyyy = cbuffer.data(fg_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxy_xyyz = cbuffer.data(fg_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxy_xyzz = cbuffer.data(fg_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xxy_xzzz = cbuffer.data(fg_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxy_yyyy = cbuffer.data(fg_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxy_yyyz = cbuffer.data(fg_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxy_yyzz = cbuffer.data(fg_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxy_yzzz = cbuffer.data(fg_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxy_zzzz = cbuffer.data(fg_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_xxz_xxxx = cbuffer.data(fg_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxz_xxxy = cbuffer.data(fg_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxz_xxxz = cbuffer.data(fg_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxz_xxyy = cbuffer.data(fg_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxz_xxyz = cbuffer.data(fg_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxz_xxzz = cbuffer.data(fg_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xxz_xyyy = cbuffer.data(fg_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xxz_xyyz = cbuffer.data(fg_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xxz_xyzz = cbuffer.data(fg_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xxz_xzzz = cbuffer.data(fg_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xxz_yyyy = cbuffer.data(fg_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xxz_yyyz = cbuffer.data(fg_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_xxz_yyzz = cbuffer.data(fg_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xxz_yzzz = cbuffer.data(fg_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xxz_zzzz = cbuffer.data(fg_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_xyy_xxxx = cbuffer.data(fg_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xyy_xxxy = cbuffer.data(fg_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xyy_xxxz = cbuffer.data(fg_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_xyy_xxyy = cbuffer.data(fg_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xyy_xxyz = cbuffer.data(fg_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_xyy_xxzz = cbuffer.data(fg_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xyy_xyyy = cbuffer.data(fg_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xyy_xyyz = cbuffer.data(fg_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xyy_xyzz = cbuffer.data(fg_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_xyy_xzzz = cbuffer.data(fg_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xyy_yyyy = cbuffer.data(fg_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xyy_yyyz = cbuffer.data(fg_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xyy_yyzz = cbuffer.data(fg_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xyy_yzzz = cbuffer.data(fg_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xyy_zzzz = cbuffer.data(fg_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_xyz_xxxx = cbuffer.data(fg_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xyz_xxxy = cbuffer.data(fg_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xyz_xxxz = cbuffer.data(fg_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xyz_xxyy = cbuffer.data(fg_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xyz_xxyz = cbuffer.data(fg_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xyz_xxzz = cbuffer.data(fg_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xyz_xyyy = cbuffer.data(fg_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xyz_xyyz = cbuffer.data(fg_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xyz_xyzz = cbuffer.data(fg_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xyz_xzzz = cbuffer.data(fg_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xyz_yyyy = cbuffer.data(fg_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xyz_yyyz = cbuffer.data(fg_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xyz_yyzz = cbuffer.data(fg_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xyz_yzzz = cbuffer.data(fg_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xyz_zzzz = cbuffer.data(fg_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xzz_xxxx = cbuffer.data(fg_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xzz_xxxy = cbuffer.data(fg_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xzz_xxxz = cbuffer.data(fg_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xzz_xxyy = cbuffer.data(fg_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xzz_xxyz = cbuffer.data(fg_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xzz_xxzz = cbuffer.data(fg_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xzz_xyyy = cbuffer.data(fg_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xzz_xyyz = cbuffer.data(fg_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xzz_xyzz = cbuffer.data(fg_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xzz_xzzz = cbuffer.data(fg_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xzz_yyyy = cbuffer.data(fg_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xzz_yyyz = cbuffer.data(fg_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xzz_yyzz = cbuffer.data(fg_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xzz_yzzz = cbuffer.data(fg_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xzz_zzzz = cbuffer.data(fg_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_yyy_xxxx = cbuffer.data(fg_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_yyy_xxxy = cbuffer.data(fg_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_yyy_xxxz = cbuffer.data(fg_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_yyy_xxyy = cbuffer.data(fg_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_yyy_xxyz = cbuffer.data(fg_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_yyy_xxzz = cbuffer.data(fg_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_yyy_xyyy = cbuffer.data(fg_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_yyy_xyyz = cbuffer.data(fg_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_yyy_xyzz = cbuffer.data(fg_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_yyy_xzzz = cbuffer.data(fg_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_yyy_yyyy = cbuffer.data(fg_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_yyy_yyyz = cbuffer.data(fg_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_yyy_yyzz = cbuffer.data(fg_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_yyy_yzzz = cbuffer.data(fg_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_yyy_zzzz = cbuffer.data(fg_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_yyz_xxxx = cbuffer.data(fg_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_yyz_xxxy = cbuffer.data(fg_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_yyz_xxxz = cbuffer.data(fg_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_yyz_xxyy = cbuffer.data(fg_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_yyz_xxyz = cbuffer.data(fg_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_yyz_xxzz = cbuffer.data(fg_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_yyz_xyyy = cbuffer.data(fg_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_yyz_xyyz = cbuffer.data(fg_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_yyz_xyzz = cbuffer.data(fg_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_yyz_xzzz = cbuffer.data(fg_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_yyz_yyyy = cbuffer.data(fg_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_yyz_yyyz = cbuffer.data(fg_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_yyz_yyzz = cbuffer.data(fg_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_yyz_yzzz = cbuffer.data(fg_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_yyz_zzzz = cbuffer.data(fg_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_yzz_xxxx = cbuffer.data(fg_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_yzz_xxxy = cbuffer.data(fg_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_yzz_xxxz = cbuffer.data(fg_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_yzz_xxyy = cbuffer.data(fg_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_yzz_xxyz = cbuffer.data(fg_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_yzz_xxzz = cbuffer.data(fg_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_yzz_xyyy = cbuffer.data(fg_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_yzz_xyyz = cbuffer.data(fg_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_yzz_xyzz = cbuffer.data(fg_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_yzz_xzzz = cbuffer.data(fg_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_yzz_yyyy = cbuffer.data(fg_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_yzz_yyyz = cbuffer.data(fg_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_yzz_yyzz = cbuffer.data(fg_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_yzz_yzzz = cbuffer.data(fg_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_yzz_zzzz = cbuffer.data(fg_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_zzz_xxxx = cbuffer.data(fg_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_zzz_xxxy = cbuffer.data(fg_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_zzz_xxxz = cbuffer.data(fg_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_zzz_xxyy = cbuffer.data(fg_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_zzz_xxyz = cbuffer.data(fg_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_zzz_xxzz = cbuffer.data(fg_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_zzz_xyyy = cbuffer.data(fg_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_zzz_xyyz = cbuffer.data(fg_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_zzz_xyzz = cbuffer.data(fg_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_zzz_xzzz = cbuffer.data(fg_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_zzz_yyyy = cbuffer.data(fg_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_zzz_yyyz = cbuffer.data(fg_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_zzz_yyzz = cbuffer.data(fg_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_zzz_yzzz = cbuffer.data(fg_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_zzz_zzzz = cbuffer.data(fg_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_z_xxx_xxxx = cbuffer.data(fg_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_xxx_xxxy = cbuffer.data(fg_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_xxx_xxxz = cbuffer.data(fg_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_xxx_xxyy = cbuffer.data(fg_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_xxx_xxyz = cbuffer.data(fg_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_xxx_xxzz = cbuffer.data(fg_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_xxx_xyyy = cbuffer.data(fg_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_xxx_xyyz = cbuffer.data(fg_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_xxx_xyzz = cbuffer.data(fg_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_xxx_xzzz = cbuffer.data(fg_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_xxx_yyyy = cbuffer.data(fg_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_xxx_yyyz = cbuffer.data(fg_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_xxx_yyzz = cbuffer.data(fg_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_xxx_yzzz = cbuffer.data(fg_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_xxx_zzzz = cbuffer.data(fg_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_xxy_xxxx = cbuffer.data(fg_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_xxy_xxxy = cbuffer.data(fg_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_xxy_xxxz = cbuffer.data(fg_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_xxy_xxyy = cbuffer.data(fg_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_xxy_xxyz = cbuffer.data(fg_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_xxy_xxzz = cbuffer.data(fg_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_xxy_xyyy = cbuffer.data(fg_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_xxy_xyyz = cbuffer.data(fg_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_xxy_xyzz = cbuffer.data(fg_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_z_xxy_xzzz = cbuffer.data(fg_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_xxy_yyyy = cbuffer.data(fg_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_xxy_yyyz = cbuffer.data(fg_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_xxy_yyzz = cbuffer.data(fg_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_xxy_yzzz = cbuffer.data(fg_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_xxy_zzzz = cbuffer.data(fg_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_z_xxz_xxxx = cbuffer.data(fg_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_xxz_xxxy = cbuffer.data(fg_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_xxz_xxxz = cbuffer.data(fg_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_xxz_xxyy = cbuffer.data(fg_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_xxz_xxyz = cbuffer.data(fg_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_xxz_xxzz = cbuffer.data(fg_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_xxz_xyyy = cbuffer.data(fg_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xxz_xyyz = cbuffer.data(fg_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xxz_xyzz = cbuffer.data(fg_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xxz_xzzz = cbuffer.data(fg_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_xxz_yyyy = cbuffer.data(fg_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xxz_yyyz = cbuffer.data(fg_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_xxz_yyzz = cbuffer.data(fg_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_xxz_yzzz = cbuffer.data(fg_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_xxz_zzzz = cbuffer.data(fg_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_xyy_xxxx = cbuffer.data(fg_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_xyy_xxxy = cbuffer.data(fg_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_xyy_xxxz = cbuffer.data(fg_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_xyy_xxyy = cbuffer.data(fg_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_xyy_xxyz = cbuffer.data(fg_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_xyy_xxzz = cbuffer.data(fg_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_xyy_xyyy = cbuffer.data(fg_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_xyy_xyyz = cbuffer.data(fg_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_xyy_xyzz = cbuffer.data(fg_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_xyy_xzzz = cbuffer.data(fg_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_xyy_yyyy = cbuffer.data(fg_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_xyy_yyyz = cbuffer.data(fg_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_xyy_yyzz = cbuffer.data(fg_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_xyy_yzzz = cbuffer.data(fg_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_xyy_zzzz = cbuffer.data(fg_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_z_xyz_xxxx = cbuffer.data(fg_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_xyz_xxxy = cbuffer.data(fg_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_xyz_xxxz = cbuffer.data(fg_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_xyz_xxyy = cbuffer.data(fg_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_xyz_xxyz = cbuffer.data(fg_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_xyz_xxzz = cbuffer.data(fg_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_xyz_xyyy = cbuffer.data(fg_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_xyz_xyyz = cbuffer.data(fg_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_xyz_xyzz = cbuffer.data(fg_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_xyz_xzzz = cbuffer.data(fg_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_xyz_yyyy = cbuffer.data(fg_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_xyz_yyyz = cbuffer.data(fg_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_xyz_yyzz = cbuffer.data(fg_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_xyz_yzzz = cbuffer.data(fg_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_xyz_zzzz = cbuffer.data(fg_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_xzz_xxxx = cbuffer.data(fg_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_xzz_xxxy = cbuffer.data(fg_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_xzz_xxxz = cbuffer.data(fg_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_z_xzz_xxyy = cbuffer.data(fg_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_z_xzz_xxyz = cbuffer.data(fg_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_z_xzz_xxzz = cbuffer.data(fg_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_z_xzz_xyyy = cbuffer.data(fg_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_z_xzz_xyyz = cbuffer.data(fg_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_z_xzz_xyzz = cbuffer.data(fg_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_z_xzz_xzzz = cbuffer.data(fg_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_z_xzz_yyyy = cbuffer.data(fg_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_z_xzz_yyyz = cbuffer.data(fg_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_z_xzz_yyzz = cbuffer.data(fg_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_z_xzz_yzzz = cbuffer.data(fg_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_z_xzz_zzzz = cbuffer.data(fg_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_z_yyy_xxxx = cbuffer.data(fg_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_z_yyy_xxxy = cbuffer.data(fg_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_z_yyy_xxxz = cbuffer.data(fg_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_z_yyy_xxyy = cbuffer.data(fg_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_z_yyy_xxyz = cbuffer.data(fg_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_z_yyy_xxzz = cbuffer.data(fg_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_z_yyy_xyyy = cbuffer.data(fg_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_z_yyy_xyyz = cbuffer.data(fg_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_z_yyy_xyzz = cbuffer.data(fg_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_z_yyy_xzzz = cbuffer.data(fg_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_z_yyy_yyyy = cbuffer.data(fg_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_z_yyy_yyyz = cbuffer.data(fg_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_z_yyy_yyzz = cbuffer.data(fg_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_z_yyy_yzzz = cbuffer.data(fg_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_z_yyy_zzzz = cbuffer.data(fg_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_z_yyz_xxxx = cbuffer.data(fg_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_z_yyz_xxxy = cbuffer.data(fg_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_z_yyz_xxxz = cbuffer.data(fg_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_z_yyz_xxyy = cbuffer.data(fg_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_z_yyz_xxyz = cbuffer.data(fg_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_z_yyz_xxzz = cbuffer.data(fg_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_z_yyz_xyyy = cbuffer.data(fg_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_z_yyz_xyyz = cbuffer.data(fg_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_z_yyz_xyzz = cbuffer.data(fg_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_z_yyz_xzzz = cbuffer.data(fg_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_z_yyz_yyyy = cbuffer.data(fg_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_z_yyz_yyyz = cbuffer.data(fg_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_z_yyz_yyzz = cbuffer.data(fg_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_z_yyz_yzzz = cbuffer.data(fg_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_z_yyz_zzzz = cbuffer.data(fg_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_z_yzz_xxxx = cbuffer.data(fg_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_yzz_xxxy = cbuffer.data(fg_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_yzz_xxxz = cbuffer.data(fg_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_yzz_xxyy = cbuffer.data(fg_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_yzz_xxyz = cbuffer.data(fg_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_yzz_xxzz = cbuffer.data(fg_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_z_yzz_xyyy = cbuffer.data(fg_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_yzz_xyyz = cbuffer.data(fg_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_yzz_xyzz = cbuffer.data(fg_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_yzz_xzzz = cbuffer.data(fg_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_z_yzz_yyyy = cbuffer.data(fg_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_yzz_yyyz = cbuffer.data(fg_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_yzz_yyzz = cbuffer.data(fg_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_yzz_yzzz = cbuffer.data(fg_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_yzz_zzzz = cbuffer.data(fg_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_zzz_xxxx = cbuffer.data(fg_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_zzz_xxxy = cbuffer.data(fg_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_zzz_xxxz = cbuffer.data(fg_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_zzz_xxyy = cbuffer.data(fg_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_zzz_xxyz = cbuffer.data(fg_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_zzz_xxzz = cbuffer.data(fg_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_zzz_xyyy = cbuffer.data(fg_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_zzz_xyyz = cbuffer.data(fg_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_zzz_xyzz = cbuffer.data(fg_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_zzz_xzzz = cbuffer.data(fg_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_zzz_yyyy = cbuffer.data(fg_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_zzz_yyyz = cbuffer.data(fg_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_zzz_yyzz = cbuffer.data(fg_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_zzz_yzzz = cbuffer.data(fg_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_zzz_zzzz = cbuffer.data(fg_geom_01_off + 449 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gfxx

            const auto gf_geom_01_off = idx_geom_01_gfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxx_xxx = cbuffer.data(gf_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_xxy = cbuffer.data(gf_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_xxz = cbuffer.data(gf_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxx_xyy = cbuffer.data(gf_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxx_xyz = cbuffer.data(gf_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxx_xzz = cbuffer.data(gf_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxx_yyy = cbuffer.data(gf_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxx_yyz = cbuffer.data(gf_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxx_yzz = cbuffer.data(gf_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxx_zzz = cbuffer.data(gf_geom_01_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_xxx, g_0_x_xxx_xxxx, g_0_x_xxx_xxxy, g_0_x_xxx_xxxz, g_0_x_xxx_xxy, g_0_x_xxx_xxyy, g_0_x_xxx_xxyz, g_0_x_xxx_xxz, g_0_x_xxx_xxzz, g_0_x_xxx_xyy, g_0_x_xxx_xyyy, g_0_x_xxx_xyyz, g_0_x_xxx_xyz, g_0_x_xxx_xyzz, g_0_x_xxx_xzz, g_0_x_xxx_xzzz, g_0_x_xxx_yyy, g_0_x_xxx_yyz, g_0_x_xxx_yzz, g_0_x_xxx_zzz, g_0_x_xxxx_xxx, g_0_x_xxxx_xxy, g_0_x_xxxx_xxz, g_0_x_xxxx_xyy, g_0_x_xxxx_xyz, g_0_x_xxxx_xzz, g_0_x_xxxx_yyy, g_0_x_xxxx_yyz, g_0_x_xxxx_yzz, g_0_x_xxxx_zzz, g_xxx_xxx, g_xxx_xxy, g_xxx_xxz, g_xxx_xyy, g_xxx_xyz, g_xxx_xzz, g_xxx_yyy, g_xxx_yyz, g_xxx_yzz, g_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxx_xxx[k] = g_xxx_xxx[k] - g_0_x_xxx_xxx[k] * ab_x + g_0_x_xxx_xxxx[k];

                g_0_x_xxxx_xxy[k] = g_xxx_xxy[k] - g_0_x_xxx_xxy[k] * ab_x + g_0_x_xxx_xxxy[k];

                g_0_x_xxxx_xxz[k] = g_xxx_xxz[k] - g_0_x_xxx_xxz[k] * ab_x + g_0_x_xxx_xxxz[k];

                g_0_x_xxxx_xyy[k] = g_xxx_xyy[k] - g_0_x_xxx_xyy[k] * ab_x + g_0_x_xxx_xxyy[k];

                g_0_x_xxxx_xyz[k] = g_xxx_xyz[k] - g_0_x_xxx_xyz[k] * ab_x + g_0_x_xxx_xxyz[k];

                g_0_x_xxxx_xzz[k] = g_xxx_xzz[k] - g_0_x_xxx_xzz[k] * ab_x + g_0_x_xxx_xxzz[k];

                g_0_x_xxxx_yyy[k] = g_xxx_yyy[k] - g_0_x_xxx_yyy[k] * ab_x + g_0_x_xxx_xyyy[k];

                g_0_x_xxxx_yyz[k] = g_xxx_yyz[k] - g_0_x_xxx_yyz[k] * ab_x + g_0_x_xxx_xyyz[k];

                g_0_x_xxxx_yzz[k] = g_xxx_yzz[k] - g_0_x_xxx_yzz[k] * ab_x + g_0_x_xxx_xyzz[k];

                g_0_x_xxxx_zzz[k] = g_xxx_zzz[k] - g_0_x_xxx_zzz[k] * ab_x + g_0_x_xxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxy_xxx = cbuffer.data(gf_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxy_xxy = cbuffer.data(gf_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxy_xxz = cbuffer.data(gf_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxy_xyy = cbuffer.data(gf_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxy_xyz = cbuffer.data(gf_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxy_xzz = cbuffer.data(gf_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxy_yyy = cbuffer.data(gf_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxy_yyz = cbuffer.data(gf_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxy_yzz = cbuffer.data(gf_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxy_zzz = cbuffer.data(gf_geom_01_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_xxx, g_0_x_xxx_xxxy, g_0_x_xxx_xxy, g_0_x_xxx_xxyy, g_0_x_xxx_xxyz, g_0_x_xxx_xxz, g_0_x_xxx_xyy, g_0_x_xxx_xyyy, g_0_x_xxx_xyyz, g_0_x_xxx_xyz, g_0_x_xxx_xyzz, g_0_x_xxx_xzz, g_0_x_xxx_yyy, g_0_x_xxx_yyyy, g_0_x_xxx_yyyz, g_0_x_xxx_yyz, g_0_x_xxx_yyzz, g_0_x_xxx_yzz, g_0_x_xxx_yzzz, g_0_x_xxx_zzz, g_0_x_xxxy_xxx, g_0_x_xxxy_xxy, g_0_x_xxxy_xxz, g_0_x_xxxy_xyy, g_0_x_xxxy_xyz, g_0_x_xxxy_xzz, g_0_x_xxxy_yyy, g_0_x_xxxy_yyz, g_0_x_xxxy_yzz, g_0_x_xxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxy_xxx[k] = -g_0_x_xxx_xxx[k] * ab_y + g_0_x_xxx_xxxy[k];

                g_0_x_xxxy_xxy[k] = -g_0_x_xxx_xxy[k] * ab_y + g_0_x_xxx_xxyy[k];

                g_0_x_xxxy_xxz[k] = -g_0_x_xxx_xxz[k] * ab_y + g_0_x_xxx_xxyz[k];

                g_0_x_xxxy_xyy[k] = -g_0_x_xxx_xyy[k] * ab_y + g_0_x_xxx_xyyy[k];

                g_0_x_xxxy_xyz[k] = -g_0_x_xxx_xyz[k] * ab_y + g_0_x_xxx_xyyz[k];

                g_0_x_xxxy_xzz[k] = -g_0_x_xxx_xzz[k] * ab_y + g_0_x_xxx_xyzz[k];

                g_0_x_xxxy_yyy[k] = -g_0_x_xxx_yyy[k] * ab_y + g_0_x_xxx_yyyy[k];

                g_0_x_xxxy_yyz[k] = -g_0_x_xxx_yyz[k] * ab_y + g_0_x_xxx_yyyz[k];

                g_0_x_xxxy_yzz[k] = -g_0_x_xxx_yzz[k] * ab_y + g_0_x_xxx_yyzz[k];

                g_0_x_xxxy_zzz[k] = -g_0_x_xxx_zzz[k] * ab_y + g_0_x_xxx_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxz_xxx = cbuffer.data(gf_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxz_xxy = cbuffer.data(gf_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxz_xxz = cbuffer.data(gf_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxz_xyy = cbuffer.data(gf_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxz_xyz = cbuffer.data(gf_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxz_xzz = cbuffer.data(gf_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxz_yyy = cbuffer.data(gf_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxz_yyz = cbuffer.data(gf_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxz_yzz = cbuffer.data(gf_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxz_zzz = cbuffer.data(gf_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_xxx, g_0_x_xxx_xxxz, g_0_x_xxx_xxy, g_0_x_xxx_xxyz, g_0_x_xxx_xxz, g_0_x_xxx_xxzz, g_0_x_xxx_xyy, g_0_x_xxx_xyyz, g_0_x_xxx_xyz, g_0_x_xxx_xyzz, g_0_x_xxx_xzz, g_0_x_xxx_xzzz, g_0_x_xxx_yyy, g_0_x_xxx_yyyz, g_0_x_xxx_yyz, g_0_x_xxx_yyzz, g_0_x_xxx_yzz, g_0_x_xxx_yzzz, g_0_x_xxx_zzz, g_0_x_xxx_zzzz, g_0_x_xxxz_xxx, g_0_x_xxxz_xxy, g_0_x_xxxz_xxz, g_0_x_xxxz_xyy, g_0_x_xxxz_xyz, g_0_x_xxxz_xzz, g_0_x_xxxz_yyy, g_0_x_xxxz_yyz, g_0_x_xxxz_yzz, g_0_x_xxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxz_xxx[k] = -g_0_x_xxx_xxx[k] * ab_z + g_0_x_xxx_xxxz[k];

                g_0_x_xxxz_xxy[k] = -g_0_x_xxx_xxy[k] * ab_z + g_0_x_xxx_xxyz[k];

                g_0_x_xxxz_xxz[k] = -g_0_x_xxx_xxz[k] * ab_z + g_0_x_xxx_xxzz[k];

                g_0_x_xxxz_xyy[k] = -g_0_x_xxx_xyy[k] * ab_z + g_0_x_xxx_xyyz[k];

                g_0_x_xxxz_xyz[k] = -g_0_x_xxx_xyz[k] * ab_z + g_0_x_xxx_xyzz[k];

                g_0_x_xxxz_xzz[k] = -g_0_x_xxx_xzz[k] * ab_z + g_0_x_xxx_xzzz[k];

                g_0_x_xxxz_yyy[k] = -g_0_x_xxx_yyy[k] * ab_z + g_0_x_xxx_yyyz[k];

                g_0_x_xxxz_yyz[k] = -g_0_x_xxx_yyz[k] * ab_z + g_0_x_xxx_yyzz[k];

                g_0_x_xxxz_yzz[k] = -g_0_x_xxx_yzz[k] * ab_z + g_0_x_xxx_yzzz[k];

                g_0_x_xxxz_zzz[k] = -g_0_x_xxx_zzz[k] * ab_z + g_0_x_xxx_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyy_xxx = cbuffer.data(gf_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxyy_xxy = cbuffer.data(gf_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxyy_xxz = cbuffer.data(gf_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxyy_xyy = cbuffer.data(gf_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxyy_xyz = cbuffer.data(gf_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxyy_xzz = cbuffer.data(gf_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxyy_yyy = cbuffer.data(gf_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxyy_yyz = cbuffer.data(gf_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxyy_yzz = cbuffer.data(gf_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxyy_zzz = cbuffer.data(gf_geom_01_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxy_xxx, g_0_x_xxy_xxxy, g_0_x_xxy_xxy, g_0_x_xxy_xxyy, g_0_x_xxy_xxyz, g_0_x_xxy_xxz, g_0_x_xxy_xyy, g_0_x_xxy_xyyy, g_0_x_xxy_xyyz, g_0_x_xxy_xyz, g_0_x_xxy_xyzz, g_0_x_xxy_xzz, g_0_x_xxy_yyy, g_0_x_xxy_yyyy, g_0_x_xxy_yyyz, g_0_x_xxy_yyz, g_0_x_xxy_yyzz, g_0_x_xxy_yzz, g_0_x_xxy_yzzz, g_0_x_xxy_zzz, g_0_x_xxyy_xxx, g_0_x_xxyy_xxy, g_0_x_xxyy_xxz, g_0_x_xxyy_xyy, g_0_x_xxyy_xyz, g_0_x_xxyy_xzz, g_0_x_xxyy_yyy, g_0_x_xxyy_yyz, g_0_x_xxyy_yzz, g_0_x_xxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyy_xxx[k] = -g_0_x_xxy_xxx[k] * ab_y + g_0_x_xxy_xxxy[k];

                g_0_x_xxyy_xxy[k] = -g_0_x_xxy_xxy[k] * ab_y + g_0_x_xxy_xxyy[k];

                g_0_x_xxyy_xxz[k] = -g_0_x_xxy_xxz[k] * ab_y + g_0_x_xxy_xxyz[k];

                g_0_x_xxyy_xyy[k] = -g_0_x_xxy_xyy[k] * ab_y + g_0_x_xxy_xyyy[k];

                g_0_x_xxyy_xyz[k] = -g_0_x_xxy_xyz[k] * ab_y + g_0_x_xxy_xyyz[k];

                g_0_x_xxyy_xzz[k] = -g_0_x_xxy_xzz[k] * ab_y + g_0_x_xxy_xyzz[k];

                g_0_x_xxyy_yyy[k] = -g_0_x_xxy_yyy[k] * ab_y + g_0_x_xxy_yyyy[k];

                g_0_x_xxyy_yyz[k] = -g_0_x_xxy_yyz[k] * ab_y + g_0_x_xxy_yyyz[k];

                g_0_x_xxyy_yzz[k] = -g_0_x_xxy_yzz[k] * ab_y + g_0_x_xxy_yyzz[k];

                g_0_x_xxyy_zzz[k] = -g_0_x_xxy_zzz[k] * ab_y + g_0_x_xxy_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyz_xxx = cbuffer.data(gf_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxyz_xxy = cbuffer.data(gf_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxyz_xxz = cbuffer.data(gf_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxyz_xyy = cbuffer.data(gf_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxyz_xyz = cbuffer.data(gf_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxyz_xzz = cbuffer.data(gf_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxyz_yyy = cbuffer.data(gf_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxyz_yyz = cbuffer.data(gf_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxyz_yzz = cbuffer.data(gf_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxyz_zzz = cbuffer.data(gf_geom_01_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyz_xxx, g_0_x_xxyz_xxy, g_0_x_xxyz_xxz, g_0_x_xxyz_xyy, g_0_x_xxyz_xyz, g_0_x_xxyz_xzz, g_0_x_xxyz_yyy, g_0_x_xxyz_yyz, g_0_x_xxyz_yzz, g_0_x_xxyz_zzz, g_0_x_xxz_xxx, g_0_x_xxz_xxxy, g_0_x_xxz_xxy, g_0_x_xxz_xxyy, g_0_x_xxz_xxyz, g_0_x_xxz_xxz, g_0_x_xxz_xyy, g_0_x_xxz_xyyy, g_0_x_xxz_xyyz, g_0_x_xxz_xyz, g_0_x_xxz_xyzz, g_0_x_xxz_xzz, g_0_x_xxz_yyy, g_0_x_xxz_yyyy, g_0_x_xxz_yyyz, g_0_x_xxz_yyz, g_0_x_xxz_yyzz, g_0_x_xxz_yzz, g_0_x_xxz_yzzz, g_0_x_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyz_xxx[k] = -g_0_x_xxz_xxx[k] * ab_y + g_0_x_xxz_xxxy[k];

                g_0_x_xxyz_xxy[k] = -g_0_x_xxz_xxy[k] * ab_y + g_0_x_xxz_xxyy[k];

                g_0_x_xxyz_xxz[k] = -g_0_x_xxz_xxz[k] * ab_y + g_0_x_xxz_xxyz[k];

                g_0_x_xxyz_xyy[k] = -g_0_x_xxz_xyy[k] * ab_y + g_0_x_xxz_xyyy[k];

                g_0_x_xxyz_xyz[k] = -g_0_x_xxz_xyz[k] * ab_y + g_0_x_xxz_xyyz[k];

                g_0_x_xxyz_xzz[k] = -g_0_x_xxz_xzz[k] * ab_y + g_0_x_xxz_xyzz[k];

                g_0_x_xxyz_yyy[k] = -g_0_x_xxz_yyy[k] * ab_y + g_0_x_xxz_yyyy[k];

                g_0_x_xxyz_yyz[k] = -g_0_x_xxz_yyz[k] * ab_y + g_0_x_xxz_yyyz[k];

                g_0_x_xxyz_yzz[k] = -g_0_x_xxz_yzz[k] * ab_y + g_0_x_xxz_yyzz[k];

                g_0_x_xxyz_zzz[k] = -g_0_x_xxz_zzz[k] * ab_y + g_0_x_xxz_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzz_xxx = cbuffer.data(gf_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxzz_xxy = cbuffer.data(gf_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxzz_xxz = cbuffer.data(gf_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxzz_xyy = cbuffer.data(gf_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxzz_xyz = cbuffer.data(gf_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxzz_xzz = cbuffer.data(gf_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxzz_yyy = cbuffer.data(gf_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxzz_yyz = cbuffer.data(gf_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxzz_yzz = cbuffer.data(gf_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxzz_zzz = cbuffer.data(gf_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxz_xxx, g_0_x_xxz_xxxz, g_0_x_xxz_xxy, g_0_x_xxz_xxyz, g_0_x_xxz_xxz, g_0_x_xxz_xxzz, g_0_x_xxz_xyy, g_0_x_xxz_xyyz, g_0_x_xxz_xyz, g_0_x_xxz_xyzz, g_0_x_xxz_xzz, g_0_x_xxz_xzzz, g_0_x_xxz_yyy, g_0_x_xxz_yyyz, g_0_x_xxz_yyz, g_0_x_xxz_yyzz, g_0_x_xxz_yzz, g_0_x_xxz_yzzz, g_0_x_xxz_zzz, g_0_x_xxz_zzzz, g_0_x_xxzz_xxx, g_0_x_xxzz_xxy, g_0_x_xxzz_xxz, g_0_x_xxzz_xyy, g_0_x_xxzz_xyz, g_0_x_xxzz_xzz, g_0_x_xxzz_yyy, g_0_x_xxzz_yyz, g_0_x_xxzz_yzz, g_0_x_xxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzz_xxx[k] = -g_0_x_xxz_xxx[k] * ab_z + g_0_x_xxz_xxxz[k];

                g_0_x_xxzz_xxy[k] = -g_0_x_xxz_xxy[k] * ab_z + g_0_x_xxz_xxyz[k];

                g_0_x_xxzz_xxz[k] = -g_0_x_xxz_xxz[k] * ab_z + g_0_x_xxz_xxzz[k];

                g_0_x_xxzz_xyy[k] = -g_0_x_xxz_xyy[k] * ab_z + g_0_x_xxz_xyyz[k];

                g_0_x_xxzz_xyz[k] = -g_0_x_xxz_xyz[k] * ab_z + g_0_x_xxz_xyzz[k];

                g_0_x_xxzz_xzz[k] = -g_0_x_xxz_xzz[k] * ab_z + g_0_x_xxz_xzzz[k];

                g_0_x_xxzz_yyy[k] = -g_0_x_xxz_yyy[k] * ab_z + g_0_x_xxz_yyyz[k];

                g_0_x_xxzz_yyz[k] = -g_0_x_xxz_yyz[k] * ab_z + g_0_x_xxz_yyzz[k];

                g_0_x_xxzz_yzz[k] = -g_0_x_xxz_yzz[k] * ab_z + g_0_x_xxz_yzzz[k];

                g_0_x_xxzz_zzz[k] = -g_0_x_xxz_zzz[k] * ab_z + g_0_x_xxz_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyy_xxx = cbuffer.data(gf_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xyyy_xxy = cbuffer.data(gf_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xyyy_xxz = cbuffer.data(gf_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xyyy_xyy = cbuffer.data(gf_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyyy_xyz = cbuffer.data(gf_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyyy_xzz = cbuffer.data(gf_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xyyy_yyy = cbuffer.data(gf_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyyy_yyz = cbuffer.data(gf_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyyy_yzz = cbuffer.data(gf_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xyyy_zzz = cbuffer.data(gf_geom_01_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyy_xxx, g_0_x_xyy_xxxy, g_0_x_xyy_xxy, g_0_x_xyy_xxyy, g_0_x_xyy_xxyz, g_0_x_xyy_xxz, g_0_x_xyy_xyy, g_0_x_xyy_xyyy, g_0_x_xyy_xyyz, g_0_x_xyy_xyz, g_0_x_xyy_xyzz, g_0_x_xyy_xzz, g_0_x_xyy_yyy, g_0_x_xyy_yyyy, g_0_x_xyy_yyyz, g_0_x_xyy_yyz, g_0_x_xyy_yyzz, g_0_x_xyy_yzz, g_0_x_xyy_yzzz, g_0_x_xyy_zzz, g_0_x_xyyy_xxx, g_0_x_xyyy_xxy, g_0_x_xyyy_xxz, g_0_x_xyyy_xyy, g_0_x_xyyy_xyz, g_0_x_xyyy_xzz, g_0_x_xyyy_yyy, g_0_x_xyyy_yyz, g_0_x_xyyy_yzz, g_0_x_xyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyy_xxx[k] = -g_0_x_xyy_xxx[k] * ab_y + g_0_x_xyy_xxxy[k];

                g_0_x_xyyy_xxy[k] = -g_0_x_xyy_xxy[k] * ab_y + g_0_x_xyy_xxyy[k];

                g_0_x_xyyy_xxz[k] = -g_0_x_xyy_xxz[k] * ab_y + g_0_x_xyy_xxyz[k];

                g_0_x_xyyy_xyy[k] = -g_0_x_xyy_xyy[k] * ab_y + g_0_x_xyy_xyyy[k];

                g_0_x_xyyy_xyz[k] = -g_0_x_xyy_xyz[k] * ab_y + g_0_x_xyy_xyyz[k];

                g_0_x_xyyy_xzz[k] = -g_0_x_xyy_xzz[k] * ab_y + g_0_x_xyy_xyzz[k];

                g_0_x_xyyy_yyy[k] = -g_0_x_xyy_yyy[k] * ab_y + g_0_x_xyy_yyyy[k];

                g_0_x_xyyy_yyz[k] = -g_0_x_xyy_yyz[k] * ab_y + g_0_x_xyy_yyyz[k];

                g_0_x_xyyy_yzz[k] = -g_0_x_xyy_yzz[k] * ab_y + g_0_x_xyy_yyzz[k];

                g_0_x_xyyy_zzz[k] = -g_0_x_xyy_zzz[k] * ab_y + g_0_x_xyy_yzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyz_xxx = cbuffer.data(gf_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyyz_xxy = cbuffer.data(gf_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xyyz_xxz = cbuffer.data(gf_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyyz_xyy = cbuffer.data(gf_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyyz_xyz = cbuffer.data(gf_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xyyz_xzz = cbuffer.data(gf_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xyyz_yyy = cbuffer.data(gf_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xyyz_yyz = cbuffer.data(gf_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xyyz_yzz = cbuffer.data(gf_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xyyz_zzz = cbuffer.data(gf_geom_01_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyz_xxx, g_0_x_xyyz_xxy, g_0_x_xyyz_xxz, g_0_x_xyyz_xyy, g_0_x_xyyz_xyz, g_0_x_xyyz_xzz, g_0_x_xyyz_yyy, g_0_x_xyyz_yyz, g_0_x_xyyz_yzz, g_0_x_xyyz_zzz, g_0_x_xyz_xxx, g_0_x_xyz_xxxy, g_0_x_xyz_xxy, g_0_x_xyz_xxyy, g_0_x_xyz_xxyz, g_0_x_xyz_xxz, g_0_x_xyz_xyy, g_0_x_xyz_xyyy, g_0_x_xyz_xyyz, g_0_x_xyz_xyz, g_0_x_xyz_xyzz, g_0_x_xyz_xzz, g_0_x_xyz_yyy, g_0_x_xyz_yyyy, g_0_x_xyz_yyyz, g_0_x_xyz_yyz, g_0_x_xyz_yyzz, g_0_x_xyz_yzz, g_0_x_xyz_yzzz, g_0_x_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyz_xxx[k] = -g_0_x_xyz_xxx[k] * ab_y + g_0_x_xyz_xxxy[k];

                g_0_x_xyyz_xxy[k] = -g_0_x_xyz_xxy[k] * ab_y + g_0_x_xyz_xxyy[k];

                g_0_x_xyyz_xxz[k] = -g_0_x_xyz_xxz[k] * ab_y + g_0_x_xyz_xxyz[k];

                g_0_x_xyyz_xyy[k] = -g_0_x_xyz_xyy[k] * ab_y + g_0_x_xyz_xyyy[k];

                g_0_x_xyyz_xyz[k] = -g_0_x_xyz_xyz[k] * ab_y + g_0_x_xyz_xyyz[k];

                g_0_x_xyyz_xzz[k] = -g_0_x_xyz_xzz[k] * ab_y + g_0_x_xyz_xyzz[k];

                g_0_x_xyyz_yyy[k] = -g_0_x_xyz_yyy[k] * ab_y + g_0_x_xyz_yyyy[k];

                g_0_x_xyyz_yyz[k] = -g_0_x_xyz_yyz[k] * ab_y + g_0_x_xyz_yyyz[k];

                g_0_x_xyyz_yzz[k] = -g_0_x_xyz_yzz[k] * ab_y + g_0_x_xyz_yyzz[k];

                g_0_x_xyyz_zzz[k] = -g_0_x_xyz_zzz[k] * ab_y + g_0_x_xyz_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzz_xxx = cbuffer.data(gf_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xyzz_xxy = cbuffer.data(gf_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xyzz_xxz = cbuffer.data(gf_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xyzz_xyy = cbuffer.data(gf_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xyzz_xyz = cbuffer.data(gf_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xyzz_xzz = cbuffer.data(gf_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xyzz_yyy = cbuffer.data(gf_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xyzz_yyz = cbuffer.data(gf_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xyzz_yzz = cbuffer.data(gf_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xyzz_zzz = cbuffer.data(gf_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzz_xxx, g_0_x_xyzz_xxy, g_0_x_xyzz_xxz, g_0_x_xyzz_xyy, g_0_x_xyzz_xyz, g_0_x_xyzz_xzz, g_0_x_xyzz_yyy, g_0_x_xyzz_yyz, g_0_x_xyzz_yzz, g_0_x_xyzz_zzz, g_0_x_xzz_xxx, g_0_x_xzz_xxxy, g_0_x_xzz_xxy, g_0_x_xzz_xxyy, g_0_x_xzz_xxyz, g_0_x_xzz_xxz, g_0_x_xzz_xyy, g_0_x_xzz_xyyy, g_0_x_xzz_xyyz, g_0_x_xzz_xyz, g_0_x_xzz_xyzz, g_0_x_xzz_xzz, g_0_x_xzz_yyy, g_0_x_xzz_yyyy, g_0_x_xzz_yyyz, g_0_x_xzz_yyz, g_0_x_xzz_yyzz, g_0_x_xzz_yzz, g_0_x_xzz_yzzz, g_0_x_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzz_xxx[k] = -g_0_x_xzz_xxx[k] * ab_y + g_0_x_xzz_xxxy[k];

                g_0_x_xyzz_xxy[k] = -g_0_x_xzz_xxy[k] * ab_y + g_0_x_xzz_xxyy[k];

                g_0_x_xyzz_xxz[k] = -g_0_x_xzz_xxz[k] * ab_y + g_0_x_xzz_xxyz[k];

                g_0_x_xyzz_xyy[k] = -g_0_x_xzz_xyy[k] * ab_y + g_0_x_xzz_xyyy[k];

                g_0_x_xyzz_xyz[k] = -g_0_x_xzz_xyz[k] * ab_y + g_0_x_xzz_xyyz[k];

                g_0_x_xyzz_xzz[k] = -g_0_x_xzz_xzz[k] * ab_y + g_0_x_xzz_xyzz[k];

                g_0_x_xyzz_yyy[k] = -g_0_x_xzz_yyy[k] * ab_y + g_0_x_xzz_yyyy[k];

                g_0_x_xyzz_yyz[k] = -g_0_x_xzz_yyz[k] * ab_y + g_0_x_xzz_yyyz[k];

                g_0_x_xyzz_yzz[k] = -g_0_x_xzz_yzz[k] * ab_y + g_0_x_xzz_yyzz[k];

                g_0_x_xyzz_zzz[k] = -g_0_x_xzz_zzz[k] * ab_y + g_0_x_xzz_yzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzz_xxx = cbuffer.data(gf_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xzzz_xxy = cbuffer.data(gf_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xzzz_xxz = cbuffer.data(gf_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xzzz_xyy = cbuffer.data(gf_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xzzz_xyz = cbuffer.data(gf_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xzzz_xzz = cbuffer.data(gf_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xzzz_yyy = cbuffer.data(gf_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xzzz_yyz = cbuffer.data(gf_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xzzz_yzz = cbuffer.data(gf_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xzzz_zzz = cbuffer.data(gf_geom_01_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzz_xxx, g_0_x_xzz_xxxz, g_0_x_xzz_xxy, g_0_x_xzz_xxyz, g_0_x_xzz_xxz, g_0_x_xzz_xxzz, g_0_x_xzz_xyy, g_0_x_xzz_xyyz, g_0_x_xzz_xyz, g_0_x_xzz_xyzz, g_0_x_xzz_xzz, g_0_x_xzz_xzzz, g_0_x_xzz_yyy, g_0_x_xzz_yyyz, g_0_x_xzz_yyz, g_0_x_xzz_yyzz, g_0_x_xzz_yzz, g_0_x_xzz_yzzz, g_0_x_xzz_zzz, g_0_x_xzz_zzzz, g_0_x_xzzz_xxx, g_0_x_xzzz_xxy, g_0_x_xzzz_xxz, g_0_x_xzzz_xyy, g_0_x_xzzz_xyz, g_0_x_xzzz_xzz, g_0_x_xzzz_yyy, g_0_x_xzzz_yyz, g_0_x_xzzz_yzz, g_0_x_xzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzz_xxx[k] = -g_0_x_xzz_xxx[k] * ab_z + g_0_x_xzz_xxxz[k];

                g_0_x_xzzz_xxy[k] = -g_0_x_xzz_xxy[k] * ab_z + g_0_x_xzz_xxyz[k];

                g_0_x_xzzz_xxz[k] = -g_0_x_xzz_xxz[k] * ab_z + g_0_x_xzz_xxzz[k];

                g_0_x_xzzz_xyy[k] = -g_0_x_xzz_xyy[k] * ab_z + g_0_x_xzz_xyyz[k];

                g_0_x_xzzz_xyz[k] = -g_0_x_xzz_xyz[k] * ab_z + g_0_x_xzz_xyzz[k];

                g_0_x_xzzz_xzz[k] = -g_0_x_xzz_xzz[k] * ab_z + g_0_x_xzz_xzzz[k];

                g_0_x_xzzz_yyy[k] = -g_0_x_xzz_yyy[k] * ab_z + g_0_x_xzz_yyyz[k];

                g_0_x_xzzz_yyz[k] = -g_0_x_xzz_yyz[k] * ab_z + g_0_x_xzz_yyzz[k];

                g_0_x_xzzz_yzz[k] = -g_0_x_xzz_yzz[k] * ab_z + g_0_x_xzz_yzzz[k];

                g_0_x_xzzz_zzz[k] = -g_0_x_xzz_zzz[k] * ab_z + g_0_x_xzz_zzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyy_xxx = cbuffer.data(gf_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yyyy_xxy = cbuffer.data(gf_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_yyyy_xxz = cbuffer.data(gf_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yyyy_xyy = cbuffer.data(gf_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yyyy_xyz = cbuffer.data(gf_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_yyyy_xzz = cbuffer.data(gf_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_yyyy_yyy = cbuffer.data(gf_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_yyyy_yyz = cbuffer.data(gf_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_yyyy_yzz = cbuffer.data(gf_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_yyyy_zzz = cbuffer.data(gf_geom_01_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyy_xxx, g_0_x_yyy_xxxy, g_0_x_yyy_xxy, g_0_x_yyy_xxyy, g_0_x_yyy_xxyz, g_0_x_yyy_xxz, g_0_x_yyy_xyy, g_0_x_yyy_xyyy, g_0_x_yyy_xyyz, g_0_x_yyy_xyz, g_0_x_yyy_xyzz, g_0_x_yyy_xzz, g_0_x_yyy_yyy, g_0_x_yyy_yyyy, g_0_x_yyy_yyyz, g_0_x_yyy_yyz, g_0_x_yyy_yyzz, g_0_x_yyy_yzz, g_0_x_yyy_yzzz, g_0_x_yyy_zzz, g_0_x_yyyy_xxx, g_0_x_yyyy_xxy, g_0_x_yyyy_xxz, g_0_x_yyyy_xyy, g_0_x_yyyy_xyz, g_0_x_yyyy_xzz, g_0_x_yyyy_yyy, g_0_x_yyyy_yyz, g_0_x_yyyy_yzz, g_0_x_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyy_xxx[k] = -g_0_x_yyy_xxx[k] * ab_y + g_0_x_yyy_xxxy[k];

                g_0_x_yyyy_xxy[k] = -g_0_x_yyy_xxy[k] * ab_y + g_0_x_yyy_xxyy[k];

                g_0_x_yyyy_xxz[k] = -g_0_x_yyy_xxz[k] * ab_y + g_0_x_yyy_xxyz[k];

                g_0_x_yyyy_xyy[k] = -g_0_x_yyy_xyy[k] * ab_y + g_0_x_yyy_xyyy[k];

                g_0_x_yyyy_xyz[k] = -g_0_x_yyy_xyz[k] * ab_y + g_0_x_yyy_xyyz[k];

                g_0_x_yyyy_xzz[k] = -g_0_x_yyy_xzz[k] * ab_y + g_0_x_yyy_xyzz[k];

                g_0_x_yyyy_yyy[k] = -g_0_x_yyy_yyy[k] * ab_y + g_0_x_yyy_yyyy[k];

                g_0_x_yyyy_yyz[k] = -g_0_x_yyy_yyz[k] * ab_y + g_0_x_yyy_yyyz[k];

                g_0_x_yyyy_yzz[k] = -g_0_x_yyy_yzz[k] * ab_y + g_0_x_yyy_yyzz[k];

                g_0_x_yyyy_zzz[k] = -g_0_x_yyy_zzz[k] * ab_y + g_0_x_yyy_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyz_xxx = cbuffer.data(gf_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_yyyz_xxy = cbuffer.data(gf_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_yyyz_xxz = cbuffer.data(gf_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_yyyz_xyy = cbuffer.data(gf_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_yyyz_xyz = cbuffer.data(gf_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_yyyz_xzz = cbuffer.data(gf_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_yyyz_yyy = cbuffer.data(gf_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_yyyz_yyz = cbuffer.data(gf_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_yyyz_yzz = cbuffer.data(gf_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_yyyz_zzz = cbuffer.data(gf_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyz_xxx, g_0_x_yyyz_xxy, g_0_x_yyyz_xxz, g_0_x_yyyz_xyy, g_0_x_yyyz_xyz, g_0_x_yyyz_xzz, g_0_x_yyyz_yyy, g_0_x_yyyz_yyz, g_0_x_yyyz_yzz, g_0_x_yyyz_zzz, g_0_x_yyz_xxx, g_0_x_yyz_xxxy, g_0_x_yyz_xxy, g_0_x_yyz_xxyy, g_0_x_yyz_xxyz, g_0_x_yyz_xxz, g_0_x_yyz_xyy, g_0_x_yyz_xyyy, g_0_x_yyz_xyyz, g_0_x_yyz_xyz, g_0_x_yyz_xyzz, g_0_x_yyz_xzz, g_0_x_yyz_yyy, g_0_x_yyz_yyyy, g_0_x_yyz_yyyz, g_0_x_yyz_yyz, g_0_x_yyz_yyzz, g_0_x_yyz_yzz, g_0_x_yyz_yzzz, g_0_x_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyz_xxx[k] = -g_0_x_yyz_xxx[k] * ab_y + g_0_x_yyz_xxxy[k];

                g_0_x_yyyz_xxy[k] = -g_0_x_yyz_xxy[k] * ab_y + g_0_x_yyz_xxyy[k];

                g_0_x_yyyz_xxz[k] = -g_0_x_yyz_xxz[k] * ab_y + g_0_x_yyz_xxyz[k];

                g_0_x_yyyz_xyy[k] = -g_0_x_yyz_xyy[k] * ab_y + g_0_x_yyz_xyyy[k];

                g_0_x_yyyz_xyz[k] = -g_0_x_yyz_xyz[k] * ab_y + g_0_x_yyz_xyyz[k];

                g_0_x_yyyz_xzz[k] = -g_0_x_yyz_xzz[k] * ab_y + g_0_x_yyz_xyzz[k];

                g_0_x_yyyz_yyy[k] = -g_0_x_yyz_yyy[k] * ab_y + g_0_x_yyz_yyyy[k];

                g_0_x_yyyz_yyz[k] = -g_0_x_yyz_yyz[k] * ab_y + g_0_x_yyz_yyyz[k];

                g_0_x_yyyz_yzz[k] = -g_0_x_yyz_yzz[k] * ab_y + g_0_x_yyz_yyzz[k];

                g_0_x_yyyz_zzz[k] = -g_0_x_yyz_zzz[k] * ab_y + g_0_x_yyz_yzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzz_xxx = cbuffer.data(gf_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_yyzz_xxy = cbuffer.data(gf_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_yyzz_xxz = cbuffer.data(gf_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_yyzz_xyy = cbuffer.data(gf_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_yyzz_xyz = cbuffer.data(gf_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_yyzz_xzz = cbuffer.data(gf_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yyzz_yyy = cbuffer.data(gf_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yyzz_yyz = cbuffer.data(gf_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yyzz_yzz = cbuffer.data(gf_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yyzz_zzz = cbuffer.data(gf_geom_01_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzz_xxx, g_0_x_yyzz_xxy, g_0_x_yyzz_xxz, g_0_x_yyzz_xyy, g_0_x_yyzz_xyz, g_0_x_yyzz_xzz, g_0_x_yyzz_yyy, g_0_x_yyzz_yyz, g_0_x_yyzz_yzz, g_0_x_yyzz_zzz, g_0_x_yzz_xxx, g_0_x_yzz_xxxy, g_0_x_yzz_xxy, g_0_x_yzz_xxyy, g_0_x_yzz_xxyz, g_0_x_yzz_xxz, g_0_x_yzz_xyy, g_0_x_yzz_xyyy, g_0_x_yzz_xyyz, g_0_x_yzz_xyz, g_0_x_yzz_xyzz, g_0_x_yzz_xzz, g_0_x_yzz_yyy, g_0_x_yzz_yyyy, g_0_x_yzz_yyyz, g_0_x_yzz_yyz, g_0_x_yzz_yyzz, g_0_x_yzz_yzz, g_0_x_yzz_yzzz, g_0_x_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzz_xxx[k] = -g_0_x_yzz_xxx[k] * ab_y + g_0_x_yzz_xxxy[k];

                g_0_x_yyzz_xxy[k] = -g_0_x_yzz_xxy[k] * ab_y + g_0_x_yzz_xxyy[k];

                g_0_x_yyzz_xxz[k] = -g_0_x_yzz_xxz[k] * ab_y + g_0_x_yzz_xxyz[k];

                g_0_x_yyzz_xyy[k] = -g_0_x_yzz_xyy[k] * ab_y + g_0_x_yzz_xyyy[k];

                g_0_x_yyzz_xyz[k] = -g_0_x_yzz_xyz[k] * ab_y + g_0_x_yzz_xyyz[k];

                g_0_x_yyzz_xzz[k] = -g_0_x_yzz_xzz[k] * ab_y + g_0_x_yzz_xyzz[k];

                g_0_x_yyzz_yyy[k] = -g_0_x_yzz_yyy[k] * ab_y + g_0_x_yzz_yyyy[k];

                g_0_x_yyzz_yyz[k] = -g_0_x_yzz_yyz[k] * ab_y + g_0_x_yzz_yyyz[k];

                g_0_x_yyzz_yzz[k] = -g_0_x_yzz_yzz[k] * ab_y + g_0_x_yzz_yyzz[k];

                g_0_x_yyzz_zzz[k] = -g_0_x_yzz_zzz[k] * ab_y + g_0_x_yzz_yzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzz_xxx = cbuffer.data(gf_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yzzz_xxy = cbuffer.data(gf_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yzzz_xxz = cbuffer.data(gf_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yzzz_xyy = cbuffer.data(gf_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yzzz_xyz = cbuffer.data(gf_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_yzzz_xzz = cbuffer.data(gf_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_yzzz_yyy = cbuffer.data(gf_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_yzzz_yyz = cbuffer.data(gf_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_yzzz_yzz = cbuffer.data(gf_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_yzzz_zzz = cbuffer.data(gf_geom_01_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzz_xxx, g_0_x_yzzz_xxy, g_0_x_yzzz_xxz, g_0_x_yzzz_xyy, g_0_x_yzzz_xyz, g_0_x_yzzz_xzz, g_0_x_yzzz_yyy, g_0_x_yzzz_yyz, g_0_x_yzzz_yzz, g_0_x_yzzz_zzz, g_0_x_zzz_xxx, g_0_x_zzz_xxxy, g_0_x_zzz_xxy, g_0_x_zzz_xxyy, g_0_x_zzz_xxyz, g_0_x_zzz_xxz, g_0_x_zzz_xyy, g_0_x_zzz_xyyy, g_0_x_zzz_xyyz, g_0_x_zzz_xyz, g_0_x_zzz_xyzz, g_0_x_zzz_xzz, g_0_x_zzz_yyy, g_0_x_zzz_yyyy, g_0_x_zzz_yyyz, g_0_x_zzz_yyz, g_0_x_zzz_yyzz, g_0_x_zzz_yzz, g_0_x_zzz_yzzz, g_0_x_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzz_xxx[k] = -g_0_x_zzz_xxx[k] * ab_y + g_0_x_zzz_xxxy[k];

                g_0_x_yzzz_xxy[k] = -g_0_x_zzz_xxy[k] * ab_y + g_0_x_zzz_xxyy[k];

                g_0_x_yzzz_xxz[k] = -g_0_x_zzz_xxz[k] * ab_y + g_0_x_zzz_xxyz[k];

                g_0_x_yzzz_xyy[k] = -g_0_x_zzz_xyy[k] * ab_y + g_0_x_zzz_xyyy[k];

                g_0_x_yzzz_xyz[k] = -g_0_x_zzz_xyz[k] * ab_y + g_0_x_zzz_xyyz[k];

                g_0_x_yzzz_xzz[k] = -g_0_x_zzz_xzz[k] * ab_y + g_0_x_zzz_xyzz[k];

                g_0_x_yzzz_yyy[k] = -g_0_x_zzz_yyy[k] * ab_y + g_0_x_zzz_yyyy[k];

                g_0_x_yzzz_yyz[k] = -g_0_x_zzz_yyz[k] * ab_y + g_0_x_zzz_yyyz[k];

                g_0_x_yzzz_yzz[k] = -g_0_x_zzz_yzz[k] * ab_y + g_0_x_zzz_yyzz[k];

                g_0_x_yzzz_zzz[k] = -g_0_x_zzz_zzz[k] * ab_y + g_0_x_zzz_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzz_xxx = cbuffer.data(gf_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_zzzz_xxy = cbuffer.data(gf_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_zzzz_xxz = cbuffer.data(gf_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_zzzz_xyy = cbuffer.data(gf_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_zzzz_xyz = cbuffer.data(gf_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_zzzz_xzz = cbuffer.data(gf_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_zzzz_yyy = cbuffer.data(gf_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_zzzz_yyz = cbuffer.data(gf_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_zzzz_yzz = cbuffer.data(gf_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_zzzz_zzz = cbuffer.data(gf_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzz_xxx, g_0_x_zzz_xxxz, g_0_x_zzz_xxy, g_0_x_zzz_xxyz, g_0_x_zzz_xxz, g_0_x_zzz_xxzz, g_0_x_zzz_xyy, g_0_x_zzz_xyyz, g_0_x_zzz_xyz, g_0_x_zzz_xyzz, g_0_x_zzz_xzz, g_0_x_zzz_xzzz, g_0_x_zzz_yyy, g_0_x_zzz_yyyz, g_0_x_zzz_yyz, g_0_x_zzz_yyzz, g_0_x_zzz_yzz, g_0_x_zzz_yzzz, g_0_x_zzz_zzz, g_0_x_zzz_zzzz, g_0_x_zzzz_xxx, g_0_x_zzzz_xxy, g_0_x_zzzz_xxz, g_0_x_zzzz_xyy, g_0_x_zzzz_xyz, g_0_x_zzzz_xzz, g_0_x_zzzz_yyy, g_0_x_zzzz_yyz, g_0_x_zzzz_yzz, g_0_x_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzz_xxx[k] = -g_0_x_zzz_xxx[k] * ab_z + g_0_x_zzz_xxxz[k];

                g_0_x_zzzz_xxy[k] = -g_0_x_zzz_xxy[k] * ab_z + g_0_x_zzz_xxyz[k];

                g_0_x_zzzz_xxz[k] = -g_0_x_zzz_xxz[k] * ab_z + g_0_x_zzz_xxzz[k];

                g_0_x_zzzz_xyy[k] = -g_0_x_zzz_xyy[k] * ab_z + g_0_x_zzz_xyyz[k];

                g_0_x_zzzz_xyz[k] = -g_0_x_zzz_xyz[k] * ab_z + g_0_x_zzz_xyzz[k];

                g_0_x_zzzz_xzz[k] = -g_0_x_zzz_xzz[k] * ab_z + g_0_x_zzz_xzzz[k];

                g_0_x_zzzz_yyy[k] = -g_0_x_zzz_yyy[k] * ab_z + g_0_x_zzz_yyyz[k];

                g_0_x_zzzz_yyz[k] = -g_0_x_zzz_yyz[k] * ab_z + g_0_x_zzz_yyzz[k];

                g_0_x_zzzz_yzz[k] = -g_0_x_zzz_yzz[k] * ab_z + g_0_x_zzz_yzzz[k];

                g_0_x_zzzz_zzz[k] = -g_0_x_zzz_zzz[k] * ab_z + g_0_x_zzz_zzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxx_xxx = cbuffer.data(gf_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xxxx_xxy = cbuffer.data(gf_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xxxx_xxz = cbuffer.data(gf_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xxxx_xyy = cbuffer.data(gf_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xxxx_xyz = cbuffer.data(gf_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xxxx_xzz = cbuffer.data(gf_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xxxx_yyy = cbuffer.data(gf_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xxxx_yyz = cbuffer.data(gf_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xxxx_yzz = cbuffer.data(gf_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xxxx_zzz = cbuffer.data(gf_geom_01_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxx_xxx, g_0_y_xxx_xxxx, g_0_y_xxx_xxxy, g_0_y_xxx_xxxz, g_0_y_xxx_xxy, g_0_y_xxx_xxyy, g_0_y_xxx_xxyz, g_0_y_xxx_xxz, g_0_y_xxx_xxzz, g_0_y_xxx_xyy, g_0_y_xxx_xyyy, g_0_y_xxx_xyyz, g_0_y_xxx_xyz, g_0_y_xxx_xyzz, g_0_y_xxx_xzz, g_0_y_xxx_xzzz, g_0_y_xxx_yyy, g_0_y_xxx_yyz, g_0_y_xxx_yzz, g_0_y_xxx_zzz, g_0_y_xxxx_xxx, g_0_y_xxxx_xxy, g_0_y_xxxx_xxz, g_0_y_xxxx_xyy, g_0_y_xxxx_xyz, g_0_y_xxxx_xzz, g_0_y_xxxx_yyy, g_0_y_xxxx_yyz, g_0_y_xxxx_yzz, g_0_y_xxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxx_xxx[k] = -g_0_y_xxx_xxx[k] * ab_x + g_0_y_xxx_xxxx[k];

                g_0_y_xxxx_xxy[k] = -g_0_y_xxx_xxy[k] * ab_x + g_0_y_xxx_xxxy[k];

                g_0_y_xxxx_xxz[k] = -g_0_y_xxx_xxz[k] * ab_x + g_0_y_xxx_xxxz[k];

                g_0_y_xxxx_xyy[k] = -g_0_y_xxx_xyy[k] * ab_x + g_0_y_xxx_xxyy[k];

                g_0_y_xxxx_xyz[k] = -g_0_y_xxx_xyz[k] * ab_x + g_0_y_xxx_xxyz[k];

                g_0_y_xxxx_xzz[k] = -g_0_y_xxx_xzz[k] * ab_x + g_0_y_xxx_xxzz[k];

                g_0_y_xxxx_yyy[k] = -g_0_y_xxx_yyy[k] * ab_x + g_0_y_xxx_xyyy[k];

                g_0_y_xxxx_yyz[k] = -g_0_y_xxx_yyz[k] * ab_x + g_0_y_xxx_xyyz[k];

                g_0_y_xxxx_yzz[k] = -g_0_y_xxx_yzz[k] * ab_x + g_0_y_xxx_xyzz[k];

                g_0_y_xxxx_zzz[k] = -g_0_y_xxx_zzz[k] * ab_x + g_0_y_xxx_xzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxy_xxx = cbuffer.data(gf_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xxxy_xxy = cbuffer.data(gf_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_xxxy_xxz = cbuffer.data(gf_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xxxy_xyy = cbuffer.data(gf_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xxxy_xyz = cbuffer.data(gf_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_xxxy_xzz = cbuffer.data(gf_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xxxy_yyy = cbuffer.data(gf_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xxxy_yyz = cbuffer.data(gf_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_xxxy_yzz = cbuffer.data(gf_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxxy_zzz = cbuffer.data(gf_geom_01_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxy_xxx, g_0_y_xxxy_xxy, g_0_y_xxxy_xxz, g_0_y_xxxy_xyy, g_0_y_xxxy_xyz, g_0_y_xxxy_xzz, g_0_y_xxxy_yyy, g_0_y_xxxy_yyz, g_0_y_xxxy_yzz, g_0_y_xxxy_zzz, g_0_y_xxy_xxx, g_0_y_xxy_xxxx, g_0_y_xxy_xxxy, g_0_y_xxy_xxxz, g_0_y_xxy_xxy, g_0_y_xxy_xxyy, g_0_y_xxy_xxyz, g_0_y_xxy_xxz, g_0_y_xxy_xxzz, g_0_y_xxy_xyy, g_0_y_xxy_xyyy, g_0_y_xxy_xyyz, g_0_y_xxy_xyz, g_0_y_xxy_xyzz, g_0_y_xxy_xzz, g_0_y_xxy_xzzz, g_0_y_xxy_yyy, g_0_y_xxy_yyz, g_0_y_xxy_yzz, g_0_y_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxy_xxx[k] = -g_0_y_xxy_xxx[k] * ab_x + g_0_y_xxy_xxxx[k];

                g_0_y_xxxy_xxy[k] = -g_0_y_xxy_xxy[k] * ab_x + g_0_y_xxy_xxxy[k];

                g_0_y_xxxy_xxz[k] = -g_0_y_xxy_xxz[k] * ab_x + g_0_y_xxy_xxxz[k];

                g_0_y_xxxy_xyy[k] = -g_0_y_xxy_xyy[k] * ab_x + g_0_y_xxy_xxyy[k];

                g_0_y_xxxy_xyz[k] = -g_0_y_xxy_xyz[k] * ab_x + g_0_y_xxy_xxyz[k];

                g_0_y_xxxy_xzz[k] = -g_0_y_xxy_xzz[k] * ab_x + g_0_y_xxy_xxzz[k];

                g_0_y_xxxy_yyy[k] = -g_0_y_xxy_yyy[k] * ab_x + g_0_y_xxy_xyyy[k];

                g_0_y_xxxy_yyz[k] = -g_0_y_xxy_yyz[k] * ab_x + g_0_y_xxy_xyyz[k];

                g_0_y_xxxy_yzz[k] = -g_0_y_xxy_yzz[k] * ab_x + g_0_y_xxy_xyzz[k];

                g_0_y_xxxy_zzz[k] = -g_0_y_xxy_zzz[k] * ab_x + g_0_y_xxy_xzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxz_xxx = cbuffer.data(gf_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxxz_xxy = cbuffer.data(gf_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxxz_xxz = cbuffer.data(gf_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxxz_xyy = cbuffer.data(gf_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xxxz_xyz = cbuffer.data(gf_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxxz_xzz = cbuffer.data(gf_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxxz_yyy = cbuffer.data(gf_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxxz_yyz = cbuffer.data(gf_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxxz_yzz = cbuffer.data(gf_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxxz_zzz = cbuffer.data(gf_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxz_xxx, g_0_y_xxxz_xxy, g_0_y_xxxz_xxz, g_0_y_xxxz_xyy, g_0_y_xxxz_xyz, g_0_y_xxxz_xzz, g_0_y_xxxz_yyy, g_0_y_xxxz_yyz, g_0_y_xxxz_yzz, g_0_y_xxxz_zzz, g_0_y_xxz_xxx, g_0_y_xxz_xxxx, g_0_y_xxz_xxxy, g_0_y_xxz_xxxz, g_0_y_xxz_xxy, g_0_y_xxz_xxyy, g_0_y_xxz_xxyz, g_0_y_xxz_xxz, g_0_y_xxz_xxzz, g_0_y_xxz_xyy, g_0_y_xxz_xyyy, g_0_y_xxz_xyyz, g_0_y_xxz_xyz, g_0_y_xxz_xyzz, g_0_y_xxz_xzz, g_0_y_xxz_xzzz, g_0_y_xxz_yyy, g_0_y_xxz_yyz, g_0_y_xxz_yzz, g_0_y_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxz_xxx[k] = -g_0_y_xxz_xxx[k] * ab_x + g_0_y_xxz_xxxx[k];

                g_0_y_xxxz_xxy[k] = -g_0_y_xxz_xxy[k] * ab_x + g_0_y_xxz_xxxy[k];

                g_0_y_xxxz_xxz[k] = -g_0_y_xxz_xxz[k] * ab_x + g_0_y_xxz_xxxz[k];

                g_0_y_xxxz_xyy[k] = -g_0_y_xxz_xyy[k] * ab_x + g_0_y_xxz_xxyy[k];

                g_0_y_xxxz_xyz[k] = -g_0_y_xxz_xyz[k] * ab_x + g_0_y_xxz_xxyz[k];

                g_0_y_xxxz_xzz[k] = -g_0_y_xxz_xzz[k] * ab_x + g_0_y_xxz_xxzz[k];

                g_0_y_xxxz_yyy[k] = -g_0_y_xxz_yyy[k] * ab_x + g_0_y_xxz_xyyy[k];

                g_0_y_xxxz_yyz[k] = -g_0_y_xxz_yyz[k] * ab_x + g_0_y_xxz_xyyz[k];

                g_0_y_xxxz_yzz[k] = -g_0_y_xxz_yzz[k] * ab_x + g_0_y_xxz_xyzz[k];

                g_0_y_xxxz_zzz[k] = -g_0_y_xxz_zzz[k] * ab_x + g_0_y_xxz_xzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyy_xxx = cbuffer.data(gf_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxyy_xxy = cbuffer.data(gf_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxyy_xxz = cbuffer.data(gf_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxyy_xyy = cbuffer.data(gf_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxyy_xyz = cbuffer.data(gf_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxyy_xzz = cbuffer.data(gf_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xxyy_yyy = cbuffer.data(gf_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xxyy_yyz = cbuffer.data(gf_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xxyy_yzz = cbuffer.data(gf_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xxyy_zzz = cbuffer.data(gf_geom_01_off + 189 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyy_xxx, g_0_y_xxyy_xxy, g_0_y_xxyy_xxz, g_0_y_xxyy_xyy, g_0_y_xxyy_xyz, g_0_y_xxyy_xzz, g_0_y_xxyy_yyy, g_0_y_xxyy_yyz, g_0_y_xxyy_yzz, g_0_y_xxyy_zzz, g_0_y_xyy_xxx, g_0_y_xyy_xxxx, g_0_y_xyy_xxxy, g_0_y_xyy_xxxz, g_0_y_xyy_xxy, g_0_y_xyy_xxyy, g_0_y_xyy_xxyz, g_0_y_xyy_xxz, g_0_y_xyy_xxzz, g_0_y_xyy_xyy, g_0_y_xyy_xyyy, g_0_y_xyy_xyyz, g_0_y_xyy_xyz, g_0_y_xyy_xyzz, g_0_y_xyy_xzz, g_0_y_xyy_xzzz, g_0_y_xyy_yyy, g_0_y_xyy_yyz, g_0_y_xyy_yzz, g_0_y_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyy_xxx[k] = -g_0_y_xyy_xxx[k] * ab_x + g_0_y_xyy_xxxx[k];

                g_0_y_xxyy_xxy[k] = -g_0_y_xyy_xxy[k] * ab_x + g_0_y_xyy_xxxy[k];

                g_0_y_xxyy_xxz[k] = -g_0_y_xyy_xxz[k] * ab_x + g_0_y_xyy_xxxz[k];

                g_0_y_xxyy_xyy[k] = -g_0_y_xyy_xyy[k] * ab_x + g_0_y_xyy_xxyy[k];

                g_0_y_xxyy_xyz[k] = -g_0_y_xyy_xyz[k] * ab_x + g_0_y_xyy_xxyz[k];

                g_0_y_xxyy_xzz[k] = -g_0_y_xyy_xzz[k] * ab_x + g_0_y_xyy_xxzz[k];

                g_0_y_xxyy_yyy[k] = -g_0_y_xyy_yyy[k] * ab_x + g_0_y_xyy_xyyy[k];

                g_0_y_xxyy_yyz[k] = -g_0_y_xyy_yyz[k] * ab_x + g_0_y_xyy_xyyz[k];

                g_0_y_xxyy_yzz[k] = -g_0_y_xyy_yzz[k] * ab_x + g_0_y_xyy_xyzz[k];

                g_0_y_xxyy_zzz[k] = -g_0_y_xyy_zzz[k] * ab_x + g_0_y_xyy_xzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyz_xxx = cbuffer.data(gf_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xxyz_xxy = cbuffer.data(gf_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_xxyz_xxz = cbuffer.data(gf_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xxyz_xyy = cbuffer.data(gf_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xxyz_xyz = cbuffer.data(gf_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_xxyz_xzz = cbuffer.data(gf_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xxyz_yyy = cbuffer.data(gf_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xxyz_yyz = cbuffer.data(gf_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_xxyz_yzz = cbuffer.data(gf_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xxyz_zzz = cbuffer.data(gf_geom_01_off + 199 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyz_xxx, g_0_y_xxyz_xxy, g_0_y_xxyz_xxz, g_0_y_xxyz_xyy, g_0_y_xxyz_xyz, g_0_y_xxyz_xzz, g_0_y_xxyz_yyy, g_0_y_xxyz_yyz, g_0_y_xxyz_yzz, g_0_y_xxyz_zzz, g_0_y_xyz_xxx, g_0_y_xyz_xxxx, g_0_y_xyz_xxxy, g_0_y_xyz_xxxz, g_0_y_xyz_xxy, g_0_y_xyz_xxyy, g_0_y_xyz_xxyz, g_0_y_xyz_xxz, g_0_y_xyz_xxzz, g_0_y_xyz_xyy, g_0_y_xyz_xyyy, g_0_y_xyz_xyyz, g_0_y_xyz_xyz, g_0_y_xyz_xyzz, g_0_y_xyz_xzz, g_0_y_xyz_xzzz, g_0_y_xyz_yyy, g_0_y_xyz_yyz, g_0_y_xyz_yzz, g_0_y_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyz_xxx[k] = -g_0_y_xyz_xxx[k] * ab_x + g_0_y_xyz_xxxx[k];

                g_0_y_xxyz_xxy[k] = -g_0_y_xyz_xxy[k] * ab_x + g_0_y_xyz_xxxy[k];

                g_0_y_xxyz_xxz[k] = -g_0_y_xyz_xxz[k] * ab_x + g_0_y_xyz_xxxz[k];

                g_0_y_xxyz_xyy[k] = -g_0_y_xyz_xyy[k] * ab_x + g_0_y_xyz_xxyy[k];

                g_0_y_xxyz_xyz[k] = -g_0_y_xyz_xyz[k] * ab_x + g_0_y_xyz_xxyz[k];

                g_0_y_xxyz_xzz[k] = -g_0_y_xyz_xzz[k] * ab_x + g_0_y_xyz_xxzz[k];

                g_0_y_xxyz_yyy[k] = -g_0_y_xyz_yyy[k] * ab_x + g_0_y_xyz_xyyy[k];

                g_0_y_xxyz_yyz[k] = -g_0_y_xyz_yyz[k] * ab_x + g_0_y_xyz_xyyz[k];

                g_0_y_xxyz_yzz[k] = -g_0_y_xyz_yzz[k] * ab_x + g_0_y_xyz_xyzz[k];

                g_0_y_xxyz_zzz[k] = -g_0_y_xyz_zzz[k] * ab_x + g_0_y_xyz_xzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzz_xxx = cbuffer.data(gf_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xxzz_xxy = cbuffer.data(gf_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xxzz_xxz = cbuffer.data(gf_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xxzz_xyy = cbuffer.data(gf_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_xxzz_xyz = cbuffer.data(gf_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xxzz_xzz = cbuffer.data(gf_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xxzz_yyy = cbuffer.data(gf_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xxzz_yyz = cbuffer.data(gf_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xxzz_yzz = cbuffer.data(gf_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xxzz_zzz = cbuffer.data(gf_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzz_xxx, g_0_y_xxzz_xxy, g_0_y_xxzz_xxz, g_0_y_xxzz_xyy, g_0_y_xxzz_xyz, g_0_y_xxzz_xzz, g_0_y_xxzz_yyy, g_0_y_xxzz_yyz, g_0_y_xxzz_yzz, g_0_y_xxzz_zzz, g_0_y_xzz_xxx, g_0_y_xzz_xxxx, g_0_y_xzz_xxxy, g_0_y_xzz_xxxz, g_0_y_xzz_xxy, g_0_y_xzz_xxyy, g_0_y_xzz_xxyz, g_0_y_xzz_xxz, g_0_y_xzz_xxzz, g_0_y_xzz_xyy, g_0_y_xzz_xyyy, g_0_y_xzz_xyyz, g_0_y_xzz_xyz, g_0_y_xzz_xyzz, g_0_y_xzz_xzz, g_0_y_xzz_xzzz, g_0_y_xzz_yyy, g_0_y_xzz_yyz, g_0_y_xzz_yzz, g_0_y_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzz_xxx[k] = -g_0_y_xzz_xxx[k] * ab_x + g_0_y_xzz_xxxx[k];

                g_0_y_xxzz_xxy[k] = -g_0_y_xzz_xxy[k] * ab_x + g_0_y_xzz_xxxy[k];

                g_0_y_xxzz_xxz[k] = -g_0_y_xzz_xxz[k] * ab_x + g_0_y_xzz_xxxz[k];

                g_0_y_xxzz_xyy[k] = -g_0_y_xzz_xyy[k] * ab_x + g_0_y_xzz_xxyy[k];

                g_0_y_xxzz_xyz[k] = -g_0_y_xzz_xyz[k] * ab_x + g_0_y_xzz_xxyz[k];

                g_0_y_xxzz_xzz[k] = -g_0_y_xzz_xzz[k] * ab_x + g_0_y_xzz_xxzz[k];

                g_0_y_xxzz_yyy[k] = -g_0_y_xzz_yyy[k] * ab_x + g_0_y_xzz_xyyy[k];

                g_0_y_xxzz_yyz[k] = -g_0_y_xzz_yyz[k] * ab_x + g_0_y_xzz_xyyz[k];

                g_0_y_xxzz_yzz[k] = -g_0_y_xzz_yzz[k] * ab_x + g_0_y_xzz_xyzz[k];

                g_0_y_xxzz_zzz[k] = -g_0_y_xzz_zzz[k] * ab_x + g_0_y_xzz_xzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyy_xxx = cbuffer.data(gf_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xyyy_xxy = cbuffer.data(gf_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xyyy_xxz = cbuffer.data(gf_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xyyy_xyy = cbuffer.data(gf_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xyyy_xyz = cbuffer.data(gf_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xyyy_xzz = cbuffer.data(gf_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xyyy_yyy = cbuffer.data(gf_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xyyy_yyz = cbuffer.data(gf_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xyyy_yzz = cbuffer.data(gf_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xyyy_zzz = cbuffer.data(gf_geom_01_off + 219 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyy_xxx, g_0_y_xyyy_xxy, g_0_y_xyyy_xxz, g_0_y_xyyy_xyy, g_0_y_xyyy_xyz, g_0_y_xyyy_xzz, g_0_y_xyyy_yyy, g_0_y_xyyy_yyz, g_0_y_xyyy_yzz, g_0_y_xyyy_zzz, g_0_y_yyy_xxx, g_0_y_yyy_xxxx, g_0_y_yyy_xxxy, g_0_y_yyy_xxxz, g_0_y_yyy_xxy, g_0_y_yyy_xxyy, g_0_y_yyy_xxyz, g_0_y_yyy_xxz, g_0_y_yyy_xxzz, g_0_y_yyy_xyy, g_0_y_yyy_xyyy, g_0_y_yyy_xyyz, g_0_y_yyy_xyz, g_0_y_yyy_xyzz, g_0_y_yyy_xzz, g_0_y_yyy_xzzz, g_0_y_yyy_yyy, g_0_y_yyy_yyz, g_0_y_yyy_yzz, g_0_y_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyy_xxx[k] = -g_0_y_yyy_xxx[k] * ab_x + g_0_y_yyy_xxxx[k];

                g_0_y_xyyy_xxy[k] = -g_0_y_yyy_xxy[k] * ab_x + g_0_y_yyy_xxxy[k];

                g_0_y_xyyy_xxz[k] = -g_0_y_yyy_xxz[k] * ab_x + g_0_y_yyy_xxxz[k];

                g_0_y_xyyy_xyy[k] = -g_0_y_yyy_xyy[k] * ab_x + g_0_y_yyy_xxyy[k];

                g_0_y_xyyy_xyz[k] = -g_0_y_yyy_xyz[k] * ab_x + g_0_y_yyy_xxyz[k];

                g_0_y_xyyy_xzz[k] = -g_0_y_yyy_xzz[k] * ab_x + g_0_y_yyy_xxzz[k];

                g_0_y_xyyy_yyy[k] = -g_0_y_yyy_yyy[k] * ab_x + g_0_y_yyy_xyyy[k];

                g_0_y_xyyy_yyz[k] = -g_0_y_yyy_yyz[k] * ab_x + g_0_y_yyy_xyyz[k];

                g_0_y_xyyy_yzz[k] = -g_0_y_yyy_yzz[k] * ab_x + g_0_y_yyy_xyzz[k];

                g_0_y_xyyy_zzz[k] = -g_0_y_yyy_zzz[k] * ab_x + g_0_y_yyy_xzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyz_xxx = cbuffer.data(gf_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xyyz_xxy = cbuffer.data(gf_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xyyz_xxz = cbuffer.data(gf_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xyyz_xyy = cbuffer.data(gf_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xyyz_xyz = cbuffer.data(gf_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xyyz_xzz = cbuffer.data(gf_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xyyz_yyy = cbuffer.data(gf_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xyyz_yyz = cbuffer.data(gf_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xyyz_yzz = cbuffer.data(gf_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xyyz_zzz = cbuffer.data(gf_geom_01_off + 229 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyz_xxx, g_0_y_xyyz_xxy, g_0_y_xyyz_xxz, g_0_y_xyyz_xyy, g_0_y_xyyz_xyz, g_0_y_xyyz_xzz, g_0_y_xyyz_yyy, g_0_y_xyyz_yyz, g_0_y_xyyz_yzz, g_0_y_xyyz_zzz, g_0_y_yyz_xxx, g_0_y_yyz_xxxx, g_0_y_yyz_xxxy, g_0_y_yyz_xxxz, g_0_y_yyz_xxy, g_0_y_yyz_xxyy, g_0_y_yyz_xxyz, g_0_y_yyz_xxz, g_0_y_yyz_xxzz, g_0_y_yyz_xyy, g_0_y_yyz_xyyy, g_0_y_yyz_xyyz, g_0_y_yyz_xyz, g_0_y_yyz_xyzz, g_0_y_yyz_xzz, g_0_y_yyz_xzzz, g_0_y_yyz_yyy, g_0_y_yyz_yyz, g_0_y_yyz_yzz, g_0_y_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyz_xxx[k] = -g_0_y_yyz_xxx[k] * ab_x + g_0_y_yyz_xxxx[k];

                g_0_y_xyyz_xxy[k] = -g_0_y_yyz_xxy[k] * ab_x + g_0_y_yyz_xxxy[k];

                g_0_y_xyyz_xxz[k] = -g_0_y_yyz_xxz[k] * ab_x + g_0_y_yyz_xxxz[k];

                g_0_y_xyyz_xyy[k] = -g_0_y_yyz_xyy[k] * ab_x + g_0_y_yyz_xxyy[k];

                g_0_y_xyyz_xyz[k] = -g_0_y_yyz_xyz[k] * ab_x + g_0_y_yyz_xxyz[k];

                g_0_y_xyyz_xzz[k] = -g_0_y_yyz_xzz[k] * ab_x + g_0_y_yyz_xxzz[k];

                g_0_y_xyyz_yyy[k] = -g_0_y_yyz_yyy[k] * ab_x + g_0_y_yyz_xyyy[k];

                g_0_y_xyyz_yyz[k] = -g_0_y_yyz_yyz[k] * ab_x + g_0_y_yyz_xyyz[k];

                g_0_y_xyyz_yzz[k] = -g_0_y_yyz_yzz[k] * ab_x + g_0_y_yyz_xyzz[k];

                g_0_y_xyyz_zzz[k] = -g_0_y_yyz_zzz[k] * ab_x + g_0_y_yyz_xzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzz_xxx = cbuffer.data(gf_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xyzz_xxy = cbuffer.data(gf_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xyzz_xxz = cbuffer.data(gf_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xyzz_xyy = cbuffer.data(gf_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xyzz_xyz = cbuffer.data(gf_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xyzz_xzz = cbuffer.data(gf_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xyzz_yyy = cbuffer.data(gf_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xyzz_yyz = cbuffer.data(gf_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xyzz_yzz = cbuffer.data(gf_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xyzz_zzz = cbuffer.data(gf_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzz_xxx, g_0_y_xyzz_xxy, g_0_y_xyzz_xxz, g_0_y_xyzz_xyy, g_0_y_xyzz_xyz, g_0_y_xyzz_xzz, g_0_y_xyzz_yyy, g_0_y_xyzz_yyz, g_0_y_xyzz_yzz, g_0_y_xyzz_zzz, g_0_y_yzz_xxx, g_0_y_yzz_xxxx, g_0_y_yzz_xxxy, g_0_y_yzz_xxxz, g_0_y_yzz_xxy, g_0_y_yzz_xxyy, g_0_y_yzz_xxyz, g_0_y_yzz_xxz, g_0_y_yzz_xxzz, g_0_y_yzz_xyy, g_0_y_yzz_xyyy, g_0_y_yzz_xyyz, g_0_y_yzz_xyz, g_0_y_yzz_xyzz, g_0_y_yzz_xzz, g_0_y_yzz_xzzz, g_0_y_yzz_yyy, g_0_y_yzz_yyz, g_0_y_yzz_yzz, g_0_y_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzz_xxx[k] = -g_0_y_yzz_xxx[k] * ab_x + g_0_y_yzz_xxxx[k];

                g_0_y_xyzz_xxy[k] = -g_0_y_yzz_xxy[k] * ab_x + g_0_y_yzz_xxxy[k];

                g_0_y_xyzz_xxz[k] = -g_0_y_yzz_xxz[k] * ab_x + g_0_y_yzz_xxxz[k];

                g_0_y_xyzz_xyy[k] = -g_0_y_yzz_xyy[k] * ab_x + g_0_y_yzz_xxyy[k];

                g_0_y_xyzz_xyz[k] = -g_0_y_yzz_xyz[k] * ab_x + g_0_y_yzz_xxyz[k];

                g_0_y_xyzz_xzz[k] = -g_0_y_yzz_xzz[k] * ab_x + g_0_y_yzz_xxzz[k];

                g_0_y_xyzz_yyy[k] = -g_0_y_yzz_yyy[k] * ab_x + g_0_y_yzz_xyyy[k];

                g_0_y_xyzz_yyz[k] = -g_0_y_yzz_yyz[k] * ab_x + g_0_y_yzz_xyyz[k];

                g_0_y_xyzz_yzz[k] = -g_0_y_yzz_yzz[k] * ab_x + g_0_y_yzz_xyzz[k];

                g_0_y_xyzz_zzz[k] = -g_0_y_yzz_zzz[k] * ab_x + g_0_y_yzz_xzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzz_xxx = cbuffer.data(gf_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xzzz_xxy = cbuffer.data(gf_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xzzz_xxz = cbuffer.data(gf_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xzzz_xyy = cbuffer.data(gf_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xzzz_xyz = cbuffer.data(gf_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xzzz_xzz = cbuffer.data(gf_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_xzzz_yyy = cbuffer.data(gf_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xzzz_yyz = cbuffer.data(gf_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xzzz_yzz = cbuffer.data(gf_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xzzz_zzz = cbuffer.data(gf_geom_01_off + 249 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzz_xxx, g_0_y_xzzz_xxy, g_0_y_xzzz_xxz, g_0_y_xzzz_xyy, g_0_y_xzzz_xyz, g_0_y_xzzz_xzz, g_0_y_xzzz_yyy, g_0_y_xzzz_yyz, g_0_y_xzzz_yzz, g_0_y_xzzz_zzz, g_0_y_zzz_xxx, g_0_y_zzz_xxxx, g_0_y_zzz_xxxy, g_0_y_zzz_xxxz, g_0_y_zzz_xxy, g_0_y_zzz_xxyy, g_0_y_zzz_xxyz, g_0_y_zzz_xxz, g_0_y_zzz_xxzz, g_0_y_zzz_xyy, g_0_y_zzz_xyyy, g_0_y_zzz_xyyz, g_0_y_zzz_xyz, g_0_y_zzz_xyzz, g_0_y_zzz_xzz, g_0_y_zzz_xzzz, g_0_y_zzz_yyy, g_0_y_zzz_yyz, g_0_y_zzz_yzz, g_0_y_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzz_xxx[k] = -g_0_y_zzz_xxx[k] * ab_x + g_0_y_zzz_xxxx[k];

                g_0_y_xzzz_xxy[k] = -g_0_y_zzz_xxy[k] * ab_x + g_0_y_zzz_xxxy[k];

                g_0_y_xzzz_xxz[k] = -g_0_y_zzz_xxz[k] * ab_x + g_0_y_zzz_xxxz[k];

                g_0_y_xzzz_xyy[k] = -g_0_y_zzz_xyy[k] * ab_x + g_0_y_zzz_xxyy[k];

                g_0_y_xzzz_xyz[k] = -g_0_y_zzz_xyz[k] * ab_x + g_0_y_zzz_xxyz[k];

                g_0_y_xzzz_xzz[k] = -g_0_y_zzz_xzz[k] * ab_x + g_0_y_zzz_xxzz[k];

                g_0_y_xzzz_yyy[k] = -g_0_y_zzz_yyy[k] * ab_x + g_0_y_zzz_xyyy[k];

                g_0_y_xzzz_yyz[k] = -g_0_y_zzz_yyz[k] * ab_x + g_0_y_zzz_xyyz[k];

                g_0_y_xzzz_yzz[k] = -g_0_y_zzz_yzz[k] * ab_x + g_0_y_zzz_xyzz[k];

                g_0_y_xzzz_zzz[k] = -g_0_y_zzz_zzz[k] * ab_x + g_0_y_zzz_xzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyy_xxx = cbuffer.data(gf_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_yyyy_xxy = cbuffer.data(gf_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_yyyy_xxz = cbuffer.data(gf_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_yyyy_xyy = cbuffer.data(gf_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_yyyy_xyz = cbuffer.data(gf_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_yyyy_xzz = cbuffer.data(gf_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_yyyy_yyy = cbuffer.data(gf_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_yyyy_yyz = cbuffer.data(gf_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_yyyy_yzz = cbuffer.data(gf_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_yyyy_zzz = cbuffer.data(gf_geom_01_off + 259 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_xxx, g_0_y_yyy_xxxy, g_0_y_yyy_xxy, g_0_y_yyy_xxyy, g_0_y_yyy_xxyz, g_0_y_yyy_xxz, g_0_y_yyy_xyy, g_0_y_yyy_xyyy, g_0_y_yyy_xyyz, g_0_y_yyy_xyz, g_0_y_yyy_xyzz, g_0_y_yyy_xzz, g_0_y_yyy_yyy, g_0_y_yyy_yyyy, g_0_y_yyy_yyyz, g_0_y_yyy_yyz, g_0_y_yyy_yyzz, g_0_y_yyy_yzz, g_0_y_yyy_yzzz, g_0_y_yyy_zzz, g_0_y_yyyy_xxx, g_0_y_yyyy_xxy, g_0_y_yyyy_xxz, g_0_y_yyyy_xyy, g_0_y_yyyy_xyz, g_0_y_yyyy_xzz, g_0_y_yyyy_yyy, g_0_y_yyyy_yyz, g_0_y_yyyy_yzz, g_0_y_yyyy_zzz, g_yyy_xxx, g_yyy_xxy, g_yyy_xxz, g_yyy_xyy, g_yyy_xyz, g_yyy_xzz, g_yyy_yyy, g_yyy_yyz, g_yyy_yzz, g_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyy_xxx[k] = g_yyy_xxx[k] - g_0_y_yyy_xxx[k] * ab_y + g_0_y_yyy_xxxy[k];

                g_0_y_yyyy_xxy[k] = g_yyy_xxy[k] - g_0_y_yyy_xxy[k] * ab_y + g_0_y_yyy_xxyy[k];

                g_0_y_yyyy_xxz[k] = g_yyy_xxz[k] - g_0_y_yyy_xxz[k] * ab_y + g_0_y_yyy_xxyz[k];

                g_0_y_yyyy_xyy[k] = g_yyy_xyy[k] - g_0_y_yyy_xyy[k] * ab_y + g_0_y_yyy_xyyy[k];

                g_0_y_yyyy_xyz[k] = g_yyy_xyz[k] - g_0_y_yyy_xyz[k] * ab_y + g_0_y_yyy_xyyz[k];

                g_0_y_yyyy_xzz[k] = g_yyy_xzz[k] - g_0_y_yyy_xzz[k] * ab_y + g_0_y_yyy_xyzz[k];

                g_0_y_yyyy_yyy[k] = g_yyy_yyy[k] - g_0_y_yyy_yyy[k] * ab_y + g_0_y_yyy_yyyy[k];

                g_0_y_yyyy_yyz[k] = g_yyy_yyz[k] - g_0_y_yyy_yyz[k] * ab_y + g_0_y_yyy_yyyz[k];

                g_0_y_yyyy_yzz[k] = g_yyy_yzz[k] - g_0_y_yyy_yzz[k] * ab_y + g_0_y_yyy_yyzz[k];

                g_0_y_yyyy_zzz[k] = g_yyy_zzz[k] - g_0_y_yyy_zzz[k] * ab_y + g_0_y_yyy_yzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyz_xxx = cbuffer.data(gf_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_yyyz_xxy = cbuffer.data(gf_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_yyyz_xxz = cbuffer.data(gf_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_yyyz_xyy = cbuffer.data(gf_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_yyyz_xyz = cbuffer.data(gf_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_yyyz_xzz = cbuffer.data(gf_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_yyyz_yyy = cbuffer.data(gf_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_yyyz_yyz = cbuffer.data(gf_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_yyyz_yzz = cbuffer.data(gf_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_yyyz_zzz = cbuffer.data(gf_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_xxx, g_0_y_yyy_xxxz, g_0_y_yyy_xxy, g_0_y_yyy_xxyz, g_0_y_yyy_xxz, g_0_y_yyy_xxzz, g_0_y_yyy_xyy, g_0_y_yyy_xyyz, g_0_y_yyy_xyz, g_0_y_yyy_xyzz, g_0_y_yyy_xzz, g_0_y_yyy_xzzz, g_0_y_yyy_yyy, g_0_y_yyy_yyyz, g_0_y_yyy_yyz, g_0_y_yyy_yyzz, g_0_y_yyy_yzz, g_0_y_yyy_yzzz, g_0_y_yyy_zzz, g_0_y_yyy_zzzz, g_0_y_yyyz_xxx, g_0_y_yyyz_xxy, g_0_y_yyyz_xxz, g_0_y_yyyz_xyy, g_0_y_yyyz_xyz, g_0_y_yyyz_xzz, g_0_y_yyyz_yyy, g_0_y_yyyz_yyz, g_0_y_yyyz_yzz, g_0_y_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyz_xxx[k] = -g_0_y_yyy_xxx[k] * ab_z + g_0_y_yyy_xxxz[k];

                g_0_y_yyyz_xxy[k] = -g_0_y_yyy_xxy[k] * ab_z + g_0_y_yyy_xxyz[k];

                g_0_y_yyyz_xxz[k] = -g_0_y_yyy_xxz[k] * ab_z + g_0_y_yyy_xxzz[k];

                g_0_y_yyyz_xyy[k] = -g_0_y_yyy_xyy[k] * ab_z + g_0_y_yyy_xyyz[k];

                g_0_y_yyyz_xyz[k] = -g_0_y_yyy_xyz[k] * ab_z + g_0_y_yyy_xyzz[k];

                g_0_y_yyyz_xzz[k] = -g_0_y_yyy_xzz[k] * ab_z + g_0_y_yyy_xzzz[k];

                g_0_y_yyyz_yyy[k] = -g_0_y_yyy_yyy[k] * ab_z + g_0_y_yyy_yyyz[k];

                g_0_y_yyyz_yyz[k] = -g_0_y_yyy_yyz[k] * ab_z + g_0_y_yyy_yyzz[k];

                g_0_y_yyyz_yzz[k] = -g_0_y_yyy_yzz[k] * ab_z + g_0_y_yyy_yzzz[k];

                g_0_y_yyyz_zzz[k] = -g_0_y_yyy_zzz[k] * ab_z + g_0_y_yyy_zzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzz_xxx = cbuffer.data(gf_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_yyzz_xxy = cbuffer.data(gf_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_yyzz_xxz = cbuffer.data(gf_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_yyzz_xyy = cbuffer.data(gf_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_yyzz_xyz = cbuffer.data(gf_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_yyzz_xzz = cbuffer.data(gf_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_yyzz_yyy = cbuffer.data(gf_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_yyzz_yyz = cbuffer.data(gf_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_yyzz_yzz = cbuffer.data(gf_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_yyzz_zzz = cbuffer.data(gf_geom_01_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyz_xxx, g_0_y_yyz_xxxz, g_0_y_yyz_xxy, g_0_y_yyz_xxyz, g_0_y_yyz_xxz, g_0_y_yyz_xxzz, g_0_y_yyz_xyy, g_0_y_yyz_xyyz, g_0_y_yyz_xyz, g_0_y_yyz_xyzz, g_0_y_yyz_xzz, g_0_y_yyz_xzzz, g_0_y_yyz_yyy, g_0_y_yyz_yyyz, g_0_y_yyz_yyz, g_0_y_yyz_yyzz, g_0_y_yyz_yzz, g_0_y_yyz_yzzz, g_0_y_yyz_zzz, g_0_y_yyz_zzzz, g_0_y_yyzz_xxx, g_0_y_yyzz_xxy, g_0_y_yyzz_xxz, g_0_y_yyzz_xyy, g_0_y_yyzz_xyz, g_0_y_yyzz_xzz, g_0_y_yyzz_yyy, g_0_y_yyzz_yyz, g_0_y_yyzz_yzz, g_0_y_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzz_xxx[k] = -g_0_y_yyz_xxx[k] * ab_z + g_0_y_yyz_xxxz[k];

                g_0_y_yyzz_xxy[k] = -g_0_y_yyz_xxy[k] * ab_z + g_0_y_yyz_xxyz[k];

                g_0_y_yyzz_xxz[k] = -g_0_y_yyz_xxz[k] * ab_z + g_0_y_yyz_xxzz[k];

                g_0_y_yyzz_xyy[k] = -g_0_y_yyz_xyy[k] * ab_z + g_0_y_yyz_xyyz[k];

                g_0_y_yyzz_xyz[k] = -g_0_y_yyz_xyz[k] * ab_z + g_0_y_yyz_xyzz[k];

                g_0_y_yyzz_xzz[k] = -g_0_y_yyz_xzz[k] * ab_z + g_0_y_yyz_xzzz[k];

                g_0_y_yyzz_yyy[k] = -g_0_y_yyz_yyy[k] * ab_z + g_0_y_yyz_yyyz[k];

                g_0_y_yyzz_yyz[k] = -g_0_y_yyz_yyz[k] * ab_z + g_0_y_yyz_yyzz[k];

                g_0_y_yyzz_yzz[k] = -g_0_y_yyz_yzz[k] * ab_z + g_0_y_yyz_yzzz[k];

                g_0_y_yyzz_zzz[k] = -g_0_y_yyz_zzz[k] * ab_z + g_0_y_yyz_zzzz[k];
            }

            /// Set up 280-290 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzz_xxx = cbuffer.data(gf_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_yzzz_xxy = cbuffer.data(gf_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_yzzz_xxz = cbuffer.data(gf_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_yzzz_xyy = cbuffer.data(gf_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_yzzz_xyz = cbuffer.data(gf_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_yzzz_xzz = cbuffer.data(gf_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_yzzz_yyy = cbuffer.data(gf_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_yzzz_yyz = cbuffer.data(gf_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_yzzz_yzz = cbuffer.data(gf_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_yzzz_zzz = cbuffer.data(gf_geom_01_off + 289 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzz_xxx, g_0_y_yzz_xxxz, g_0_y_yzz_xxy, g_0_y_yzz_xxyz, g_0_y_yzz_xxz, g_0_y_yzz_xxzz, g_0_y_yzz_xyy, g_0_y_yzz_xyyz, g_0_y_yzz_xyz, g_0_y_yzz_xyzz, g_0_y_yzz_xzz, g_0_y_yzz_xzzz, g_0_y_yzz_yyy, g_0_y_yzz_yyyz, g_0_y_yzz_yyz, g_0_y_yzz_yyzz, g_0_y_yzz_yzz, g_0_y_yzz_yzzz, g_0_y_yzz_zzz, g_0_y_yzz_zzzz, g_0_y_yzzz_xxx, g_0_y_yzzz_xxy, g_0_y_yzzz_xxz, g_0_y_yzzz_xyy, g_0_y_yzzz_xyz, g_0_y_yzzz_xzz, g_0_y_yzzz_yyy, g_0_y_yzzz_yyz, g_0_y_yzzz_yzz, g_0_y_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzz_xxx[k] = -g_0_y_yzz_xxx[k] * ab_z + g_0_y_yzz_xxxz[k];

                g_0_y_yzzz_xxy[k] = -g_0_y_yzz_xxy[k] * ab_z + g_0_y_yzz_xxyz[k];

                g_0_y_yzzz_xxz[k] = -g_0_y_yzz_xxz[k] * ab_z + g_0_y_yzz_xxzz[k];

                g_0_y_yzzz_xyy[k] = -g_0_y_yzz_xyy[k] * ab_z + g_0_y_yzz_xyyz[k];

                g_0_y_yzzz_xyz[k] = -g_0_y_yzz_xyz[k] * ab_z + g_0_y_yzz_xyzz[k];

                g_0_y_yzzz_xzz[k] = -g_0_y_yzz_xzz[k] * ab_z + g_0_y_yzz_xzzz[k];

                g_0_y_yzzz_yyy[k] = -g_0_y_yzz_yyy[k] * ab_z + g_0_y_yzz_yyyz[k];

                g_0_y_yzzz_yyz[k] = -g_0_y_yzz_yyz[k] * ab_z + g_0_y_yzz_yyzz[k];

                g_0_y_yzzz_yzz[k] = -g_0_y_yzz_yzz[k] * ab_z + g_0_y_yzz_yzzz[k];

                g_0_y_yzzz_zzz[k] = -g_0_y_yzz_zzz[k] * ab_z + g_0_y_yzz_zzzz[k];
            }

            /// Set up 290-300 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzz_xxx = cbuffer.data(gf_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_zzzz_xxy = cbuffer.data(gf_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_zzzz_xxz = cbuffer.data(gf_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_zzzz_xyy = cbuffer.data(gf_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_zzzz_xyz = cbuffer.data(gf_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_zzzz_xzz = cbuffer.data(gf_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_zzzz_yyy = cbuffer.data(gf_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_zzzz_yyz = cbuffer.data(gf_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_zzzz_yzz = cbuffer.data(gf_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_zzzz_zzz = cbuffer.data(gf_geom_01_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzz_xxx, g_0_y_zzz_xxxz, g_0_y_zzz_xxy, g_0_y_zzz_xxyz, g_0_y_zzz_xxz, g_0_y_zzz_xxzz, g_0_y_zzz_xyy, g_0_y_zzz_xyyz, g_0_y_zzz_xyz, g_0_y_zzz_xyzz, g_0_y_zzz_xzz, g_0_y_zzz_xzzz, g_0_y_zzz_yyy, g_0_y_zzz_yyyz, g_0_y_zzz_yyz, g_0_y_zzz_yyzz, g_0_y_zzz_yzz, g_0_y_zzz_yzzz, g_0_y_zzz_zzz, g_0_y_zzz_zzzz, g_0_y_zzzz_xxx, g_0_y_zzzz_xxy, g_0_y_zzzz_xxz, g_0_y_zzzz_xyy, g_0_y_zzzz_xyz, g_0_y_zzzz_xzz, g_0_y_zzzz_yyy, g_0_y_zzzz_yyz, g_0_y_zzzz_yzz, g_0_y_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzz_xxx[k] = -g_0_y_zzz_xxx[k] * ab_z + g_0_y_zzz_xxxz[k];

                g_0_y_zzzz_xxy[k] = -g_0_y_zzz_xxy[k] * ab_z + g_0_y_zzz_xxyz[k];

                g_0_y_zzzz_xxz[k] = -g_0_y_zzz_xxz[k] * ab_z + g_0_y_zzz_xxzz[k];

                g_0_y_zzzz_xyy[k] = -g_0_y_zzz_xyy[k] * ab_z + g_0_y_zzz_xyyz[k];

                g_0_y_zzzz_xyz[k] = -g_0_y_zzz_xyz[k] * ab_z + g_0_y_zzz_xyzz[k];

                g_0_y_zzzz_xzz[k] = -g_0_y_zzz_xzz[k] * ab_z + g_0_y_zzz_xzzz[k];

                g_0_y_zzzz_yyy[k] = -g_0_y_zzz_yyy[k] * ab_z + g_0_y_zzz_yyyz[k];

                g_0_y_zzzz_yyz[k] = -g_0_y_zzz_yyz[k] * ab_z + g_0_y_zzz_yyzz[k];

                g_0_y_zzzz_yzz[k] = -g_0_y_zzz_yzz[k] * ab_z + g_0_y_zzz_yzzz[k];

                g_0_y_zzzz_zzz[k] = -g_0_y_zzz_zzz[k] * ab_z + g_0_y_zzz_zzzz[k];
            }

            /// Set up 300-310 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxx_xxx = cbuffer.data(gf_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_xxxx_xxy = cbuffer.data(gf_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_xxxx_xxz = cbuffer.data(gf_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_xxxx_xyy = cbuffer.data(gf_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_xxxx_xyz = cbuffer.data(gf_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_xxxx_xzz = cbuffer.data(gf_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_xxxx_yyy = cbuffer.data(gf_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_xxxx_yyz = cbuffer.data(gf_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_xxxx_yzz = cbuffer.data(gf_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_xxxx_zzz = cbuffer.data(gf_geom_01_off + 309 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxx_xxx, g_0_z_xxx_xxxx, g_0_z_xxx_xxxy, g_0_z_xxx_xxxz, g_0_z_xxx_xxy, g_0_z_xxx_xxyy, g_0_z_xxx_xxyz, g_0_z_xxx_xxz, g_0_z_xxx_xxzz, g_0_z_xxx_xyy, g_0_z_xxx_xyyy, g_0_z_xxx_xyyz, g_0_z_xxx_xyz, g_0_z_xxx_xyzz, g_0_z_xxx_xzz, g_0_z_xxx_xzzz, g_0_z_xxx_yyy, g_0_z_xxx_yyz, g_0_z_xxx_yzz, g_0_z_xxx_zzz, g_0_z_xxxx_xxx, g_0_z_xxxx_xxy, g_0_z_xxxx_xxz, g_0_z_xxxx_xyy, g_0_z_xxxx_xyz, g_0_z_xxxx_xzz, g_0_z_xxxx_yyy, g_0_z_xxxx_yyz, g_0_z_xxxx_yzz, g_0_z_xxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxx_xxx[k] = -g_0_z_xxx_xxx[k] * ab_x + g_0_z_xxx_xxxx[k];

                g_0_z_xxxx_xxy[k] = -g_0_z_xxx_xxy[k] * ab_x + g_0_z_xxx_xxxy[k];

                g_0_z_xxxx_xxz[k] = -g_0_z_xxx_xxz[k] * ab_x + g_0_z_xxx_xxxz[k];

                g_0_z_xxxx_xyy[k] = -g_0_z_xxx_xyy[k] * ab_x + g_0_z_xxx_xxyy[k];

                g_0_z_xxxx_xyz[k] = -g_0_z_xxx_xyz[k] * ab_x + g_0_z_xxx_xxyz[k];

                g_0_z_xxxx_xzz[k] = -g_0_z_xxx_xzz[k] * ab_x + g_0_z_xxx_xxzz[k];

                g_0_z_xxxx_yyy[k] = -g_0_z_xxx_yyy[k] * ab_x + g_0_z_xxx_xyyy[k];

                g_0_z_xxxx_yyz[k] = -g_0_z_xxx_yyz[k] * ab_x + g_0_z_xxx_xyyz[k];

                g_0_z_xxxx_yzz[k] = -g_0_z_xxx_yzz[k] * ab_x + g_0_z_xxx_xyzz[k];

                g_0_z_xxxx_zzz[k] = -g_0_z_xxx_zzz[k] * ab_x + g_0_z_xxx_xzzz[k];
            }

            /// Set up 310-320 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxy_xxx = cbuffer.data(gf_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_xxxy_xxy = cbuffer.data(gf_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_xxxy_xxz = cbuffer.data(gf_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_xxxy_xyy = cbuffer.data(gf_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_xxxy_xyz = cbuffer.data(gf_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_xxxy_xzz = cbuffer.data(gf_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_xxxy_yyy = cbuffer.data(gf_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_xxxy_yyz = cbuffer.data(gf_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_xxxy_yzz = cbuffer.data(gf_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_xxxy_zzz = cbuffer.data(gf_geom_01_off + 319 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxy_xxx, g_0_z_xxxy_xxy, g_0_z_xxxy_xxz, g_0_z_xxxy_xyy, g_0_z_xxxy_xyz, g_0_z_xxxy_xzz, g_0_z_xxxy_yyy, g_0_z_xxxy_yyz, g_0_z_xxxy_yzz, g_0_z_xxxy_zzz, g_0_z_xxy_xxx, g_0_z_xxy_xxxx, g_0_z_xxy_xxxy, g_0_z_xxy_xxxz, g_0_z_xxy_xxy, g_0_z_xxy_xxyy, g_0_z_xxy_xxyz, g_0_z_xxy_xxz, g_0_z_xxy_xxzz, g_0_z_xxy_xyy, g_0_z_xxy_xyyy, g_0_z_xxy_xyyz, g_0_z_xxy_xyz, g_0_z_xxy_xyzz, g_0_z_xxy_xzz, g_0_z_xxy_xzzz, g_0_z_xxy_yyy, g_0_z_xxy_yyz, g_0_z_xxy_yzz, g_0_z_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxy_xxx[k] = -g_0_z_xxy_xxx[k] * ab_x + g_0_z_xxy_xxxx[k];

                g_0_z_xxxy_xxy[k] = -g_0_z_xxy_xxy[k] * ab_x + g_0_z_xxy_xxxy[k];

                g_0_z_xxxy_xxz[k] = -g_0_z_xxy_xxz[k] * ab_x + g_0_z_xxy_xxxz[k];

                g_0_z_xxxy_xyy[k] = -g_0_z_xxy_xyy[k] * ab_x + g_0_z_xxy_xxyy[k];

                g_0_z_xxxy_xyz[k] = -g_0_z_xxy_xyz[k] * ab_x + g_0_z_xxy_xxyz[k];

                g_0_z_xxxy_xzz[k] = -g_0_z_xxy_xzz[k] * ab_x + g_0_z_xxy_xxzz[k];

                g_0_z_xxxy_yyy[k] = -g_0_z_xxy_yyy[k] * ab_x + g_0_z_xxy_xyyy[k];

                g_0_z_xxxy_yyz[k] = -g_0_z_xxy_yyz[k] * ab_x + g_0_z_xxy_xyyz[k];

                g_0_z_xxxy_yzz[k] = -g_0_z_xxy_yzz[k] * ab_x + g_0_z_xxy_xyzz[k];

                g_0_z_xxxy_zzz[k] = -g_0_z_xxy_zzz[k] * ab_x + g_0_z_xxy_xzzz[k];
            }

            /// Set up 320-330 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxz_xxx = cbuffer.data(gf_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_xxxz_xxy = cbuffer.data(gf_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_xxxz_xxz = cbuffer.data(gf_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_xxxz_xyy = cbuffer.data(gf_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_z_xxxz_xyz = cbuffer.data(gf_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_xxxz_xzz = cbuffer.data(gf_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_xxxz_yyy = cbuffer.data(gf_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_xxxz_yyz = cbuffer.data(gf_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_xxxz_yzz = cbuffer.data(gf_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_xxxz_zzz = cbuffer.data(gf_geom_01_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxz_xxx, g_0_z_xxxz_xxy, g_0_z_xxxz_xxz, g_0_z_xxxz_xyy, g_0_z_xxxz_xyz, g_0_z_xxxz_xzz, g_0_z_xxxz_yyy, g_0_z_xxxz_yyz, g_0_z_xxxz_yzz, g_0_z_xxxz_zzz, g_0_z_xxz_xxx, g_0_z_xxz_xxxx, g_0_z_xxz_xxxy, g_0_z_xxz_xxxz, g_0_z_xxz_xxy, g_0_z_xxz_xxyy, g_0_z_xxz_xxyz, g_0_z_xxz_xxz, g_0_z_xxz_xxzz, g_0_z_xxz_xyy, g_0_z_xxz_xyyy, g_0_z_xxz_xyyz, g_0_z_xxz_xyz, g_0_z_xxz_xyzz, g_0_z_xxz_xzz, g_0_z_xxz_xzzz, g_0_z_xxz_yyy, g_0_z_xxz_yyz, g_0_z_xxz_yzz, g_0_z_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxz_xxx[k] = -g_0_z_xxz_xxx[k] * ab_x + g_0_z_xxz_xxxx[k];

                g_0_z_xxxz_xxy[k] = -g_0_z_xxz_xxy[k] * ab_x + g_0_z_xxz_xxxy[k];

                g_0_z_xxxz_xxz[k] = -g_0_z_xxz_xxz[k] * ab_x + g_0_z_xxz_xxxz[k];

                g_0_z_xxxz_xyy[k] = -g_0_z_xxz_xyy[k] * ab_x + g_0_z_xxz_xxyy[k];

                g_0_z_xxxz_xyz[k] = -g_0_z_xxz_xyz[k] * ab_x + g_0_z_xxz_xxyz[k];

                g_0_z_xxxz_xzz[k] = -g_0_z_xxz_xzz[k] * ab_x + g_0_z_xxz_xxzz[k];

                g_0_z_xxxz_yyy[k] = -g_0_z_xxz_yyy[k] * ab_x + g_0_z_xxz_xyyy[k];

                g_0_z_xxxz_yyz[k] = -g_0_z_xxz_yyz[k] * ab_x + g_0_z_xxz_xyyz[k];

                g_0_z_xxxz_yzz[k] = -g_0_z_xxz_yzz[k] * ab_x + g_0_z_xxz_xyzz[k];

                g_0_z_xxxz_zzz[k] = -g_0_z_xxz_zzz[k] * ab_x + g_0_z_xxz_xzzz[k];
            }

            /// Set up 330-340 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyy_xxx = cbuffer.data(gf_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_xxyy_xxy = cbuffer.data(gf_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_xxyy_xxz = cbuffer.data(gf_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_xxyy_xyy = cbuffer.data(gf_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_xxyy_xyz = cbuffer.data(gf_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_xxyy_xzz = cbuffer.data(gf_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_xxyy_yyy = cbuffer.data(gf_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xxyy_yyz = cbuffer.data(gf_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xxyy_yzz = cbuffer.data(gf_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xxyy_zzz = cbuffer.data(gf_geom_01_off + 339 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyy_xxx, g_0_z_xxyy_xxy, g_0_z_xxyy_xxz, g_0_z_xxyy_xyy, g_0_z_xxyy_xyz, g_0_z_xxyy_xzz, g_0_z_xxyy_yyy, g_0_z_xxyy_yyz, g_0_z_xxyy_yzz, g_0_z_xxyy_zzz, g_0_z_xyy_xxx, g_0_z_xyy_xxxx, g_0_z_xyy_xxxy, g_0_z_xyy_xxxz, g_0_z_xyy_xxy, g_0_z_xyy_xxyy, g_0_z_xyy_xxyz, g_0_z_xyy_xxz, g_0_z_xyy_xxzz, g_0_z_xyy_xyy, g_0_z_xyy_xyyy, g_0_z_xyy_xyyz, g_0_z_xyy_xyz, g_0_z_xyy_xyzz, g_0_z_xyy_xzz, g_0_z_xyy_xzzz, g_0_z_xyy_yyy, g_0_z_xyy_yyz, g_0_z_xyy_yzz, g_0_z_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyy_xxx[k] = -g_0_z_xyy_xxx[k] * ab_x + g_0_z_xyy_xxxx[k];

                g_0_z_xxyy_xxy[k] = -g_0_z_xyy_xxy[k] * ab_x + g_0_z_xyy_xxxy[k];

                g_0_z_xxyy_xxz[k] = -g_0_z_xyy_xxz[k] * ab_x + g_0_z_xyy_xxxz[k];

                g_0_z_xxyy_xyy[k] = -g_0_z_xyy_xyy[k] * ab_x + g_0_z_xyy_xxyy[k];

                g_0_z_xxyy_xyz[k] = -g_0_z_xyy_xyz[k] * ab_x + g_0_z_xyy_xxyz[k];

                g_0_z_xxyy_xzz[k] = -g_0_z_xyy_xzz[k] * ab_x + g_0_z_xyy_xxzz[k];

                g_0_z_xxyy_yyy[k] = -g_0_z_xyy_yyy[k] * ab_x + g_0_z_xyy_xyyy[k];

                g_0_z_xxyy_yyz[k] = -g_0_z_xyy_yyz[k] * ab_x + g_0_z_xyy_xyyz[k];

                g_0_z_xxyy_yzz[k] = -g_0_z_xyy_yzz[k] * ab_x + g_0_z_xyy_xyzz[k];

                g_0_z_xxyy_zzz[k] = -g_0_z_xyy_zzz[k] * ab_x + g_0_z_xyy_xzzz[k];
            }

            /// Set up 340-350 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyz_xxx = cbuffer.data(gf_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xxyz_xxy = cbuffer.data(gf_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_xxyz_xxz = cbuffer.data(gf_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_xxyz_xyy = cbuffer.data(gf_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_xxyz_xyz = cbuffer.data(gf_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_xxyz_xzz = cbuffer.data(gf_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_xxyz_yyy = cbuffer.data(gf_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_xxyz_yyz = cbuffer.data(gf_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_xxyz_yzz = cbuffer.data(gf_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_xxyz_zzz = cbuffer.data(gf_geom_01_off + 349 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyz_xxx, g_0_z_xxyz_xxy, g_0_z_xxyz_xxz, g_0_z_xxyz_xyy, g_0_z_xxyz_xyz, g_0_z_xxyz_xzz, g_0_z_xxyz_yyy, g_0_z_xxyz_yyz, g_0_z_xxyz_yzz, g_0_z_xxyz_zzz, g_0_z_xyz_xxx, g_0_z_xyz_xxxx, g_0_z_xyz_xxxy, g_0_z_xyz_xxxz, g_0_z_xyz_xxy, g_0_z_xyz_xxyy, g_0_z_xyz_xxyz, g_0_z_xyz_xxz, g_0_z_xyz_xxzz, g_0_z_xyz_xyy, g_0_z_xyz_xyyy, g_0_z_xyz_xyyz, g_0_z_xyz_xyz, g_0_z_xyz_xyzz, g_0_z_xyz_xzz, g_0_z_xyz_xzzz, g_0_z_xyz_yyy, g_0_z_xyz_yyz, g_0_z_xyz_yzz, g_0_z_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyz_xxx[k] = -g_0_z_xyz_xxx[k] * ab_x + g_0_z_xyz_xxxx[k];

                g_0_z_xxyz_xxy[k] = -g_0_z_xyz_xxy[k] * ab_x + g_0_z_xyz_xxxy[k];

                g_0_z_xxyz_xxz[k] = -g_0_z_xyz_xxz[k] * ab_x + g_0_z_xyz_xxxz[k];

                g_0_z_xxyz_xyy[k] = -g_0_z_xyz_xyy[k] * ab_x + g_0_z_xyz_xxyy[k];

                g_0_z_xxyz_xyz[k] = -g_0_z_xyz_xyz[k] * ab_x + g_0_z_xyz_xxyz[k];

                g_0_z_xxyz_xzz[k] = -g_0_z_xyz_xzz[k] * ab_x + g_0_z_xyz_xxzz[k];

                g_0_z_xxyz_yyy[k] = -g_0_z_xyz_yyy[k] * ab_x + g_0_z_xyz_xyyy[k];

                g_0_z_xxyz_yyz[k] = -g_0_z_xyz_yyz[k] * ab_x + g_0_z_xyz_xyyz[k];

                g_0_z_xxyz_yzz[k] = -g_0_z_xyz_yzz[k] * ab_x + g_0_z_xyz_xyzz[k];

                g_0_z_xxyz_zzz[k] = -g_0_z_xyz_zzz[k] * ab_x + g_0_z_xyz_xzzz[k];
            }

            /// Set up 350-360 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzz_xxx = cbuffer.data(gf_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_xxzz_xxy = cbuffer.data(gf_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_xxzz_xxz = cbuffer.data(gf_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_xxzz_xyy = cbuffer.data(gf_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_xxzz_xyz = cbuffer.data(gf_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_xxzz_xzz = cbuffer.data(gf_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_xxzz_yyy = cbuffer.data(gf_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_xxzz_yyz = cbuffer.data(gf_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_xxzz_yzz = cbuffer.data(gf_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_xxzz_zzz = cbuffer.data(gf_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzz_xxx, g_0_z_xxzz_xxy, g_0_z_xxzz_xxz, g_0_z_xxzz_xyy, g_0_z_xxzz_xyz, g_0_z_xxzz_xzz, g_0_z_xxzz_yyy, g_0_z_xxzz_yyz, g_0_z_xxzz_yzz, g_0_z_xxzz_zzz, g_0_z_xzz_xxx, g_0_z_xzz_xxxx, g_0_z_xzz_xxxy, g_0_z_xzz_xxxz, g_0_z_xzz_xxy, g_0_z_xzz_xxyy, g_0_z_xzz_xxyz, g_0_z_xzz_xxz, g_0_z_xzz_xxzz, g_0_z_xzz_xyy, g_0_z_xzz_xyyy, g_0_z_xzz_xyyz, g_0_z_xzz_xyz, g_0_z_xzz_xyzz, g_0_z_xzz_xzz, g_0_z_xzz_xzzz, g_0_z_xzz_yyy, g_0_z_xzz_yyz, g_0_z_xzz_yzz, g_0_z_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzz_xxx[k] = -g_0_z_xzz_xxx[k] * ab_x + g_0_z_xzz_xxxx[k];

                g_0_z_xxzz_xxy[k] = -g_0_z_xzz_xxy[k] * ab_x + g_0_z_xzz_xxxy[k];

                g_0_z_xxzz_xxz[k] = -g_0_z_xzz_xxz[k] * ab_x + g_0_z_xzz_xxxz[k];

                g_0_z_xxzz_xyy[k] = -g_0_z_xzz_xyy[k] * ab_x + g_0_z_xzz_xxyy[k];

                g_0_z_xxzz_xyz[k] = -g_0_z_xzz_xyz[k] * ab_x + g_0_z_xzz_xxyz[k];

                g_0_z_xxzz_xzz[k] = -g_0_z_xzz_xzz[k] * ab_x + g_0_z_xzz_xxzz[k];

                g_0_z_xxzz_yyy[k] = -g_0_z_xzz_yyy[k] * ab_x + g_0_z_xzz_xyyy[k];

                g_0_z_xxzz_yyz[k] = -g_0_z_xzz_yyz[k] * ab_x + g_0_z_xzz_xyyz[k];

                g_0_z_xxzz_yzz[k] = -g_0_z_xzz_yzz[k] * ab_x + g_0_z_xzz_xyzz[k];

                g_0_z_xxzz_zzz[k] = -g_0_z_xzz_zzz[k] * ab_x + g_0_z_xzz_xzzz[k];
            }

            /// Set up 360-370 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyy_xxx = cbuffer.data(gf_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_xyyy_xxy = cbuffer.data(gf_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_xyyy_xxz = cbuffer.data(gf_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_xyyy_xyy = cbuffer.data(gf_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_xyyy_xyz = cbuffer.data(gf_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_xyyy_xzz = cbuffer.data(gf_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_xyyy_yyy = cbuffer.data(gf_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_xyyy_yyz = cbuffer.data(gf_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_xyyy_yzz = cbuffer.data(gf_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_xyyy_zzz = cbuffer.data(gf_geom_01_off + 369 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyy_xxx, g_0_z_xyyy_xxy, g_0_z_xyyy_xxz, g_0_z_xyyy_xyy, g_0_z_xyyy_xyz, g_0_z_xyyy_xzz, g_0_z_xyyy_yyy, g_0_z_xyyy_yyz, g_0_z_xyyy_yzz, g_0_z_xyyy_zzz, g_0_z_yyy_xxx, g_0_z_yyy_xxxx, g_0_z_yyy_xxxy, g_0_z_yyy_xxxz, g_0_z_yyy_xxy, g_0_z_yyy_xxyy, g_0_z_yyy_xxyz, g_0_z_yyy_xxz, g_0_z_yyy_xxzz, g_0_z_yyy_xyy, g_0_z_yyy_xyyy, g_0_z_yyy_xyyz, g_0_z_yyy_xyz, g_0_z_yyy_xyzz, g_0_z_yyy_xzz, g_0_z_yyy_xzzz, g_0_z_yyy_yyy, g_0_z_yyy_yyz, g_0_z_yyy_yzz, g_0_z_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyy_xxx[k] = -g_0_z_yyy_xxx[k] * ab_x + g_0_z_yyy_xxxx[k];

                g_0_z_xyyy_xxy[k] = -g_0_z_yyy_xxy[k] * ab_x + g_0_z_yyy_xxxy[k];

                g_0_z_xyyy_xxz[k] = -g_0_z_yyy_xxz[k] * ab_x + g_0_z_yyy_xxxz[k];

                g_0_z_xyyy_xyy[k] = -g_0_z_yyy_xyy[k] * ab_x + g_0_z_yyy_xxyy[k];

                g_0_z_xyyy_xyz[k] = -g_0_z_yyy_xyz[k] * ab_x + g_0_z_yyy_xxyz[k];

                g_0_z_xyyy_xzz[k] = -g_0_z_yyy_xzz[k] * ab_x + g_0_z_yyy_xxzz[k];

                g_0_z_xyyy_yyy[k] = -g_0_z_yyy_yyy[k] * ab_x + g_0_z_yyy_xyyy[k];

                g_0_z_xyyy_yyz[k] = -g_0_z_yyy_yyz[k] * ab_x + g_0_z_yyy_xyyz[k];

                g_0_z_xyyy_yzz[k] = -g_0_z_yyy_yzz[k] * ab_x + g_0_z_yyy_xyzz[k];

                g_0_z_xyyy_zzz[k] = -g_0_z_yyy_zzz[k] * ab_x + g_0_z_yyy_xzzz[k];
            }

            /// Set up 370-380 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyz_xxx = cbuffer.data(gf_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_xyyz_xxy = cbuffer.data(gf_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_xyyz_xxz = cbuffer.data(gf_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_xyyz_xyy = cbuffer.data(gf_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_xyyz_xyz = cbuffer.data(gf_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_xyyz_xzz = cbuffer.data(gf_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_xyyz_yyy = cbuffer.data(gf_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_xyyz_yyz = cbuffer.data(gf_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_z_xyyz_yzz = cbuffer.data(gf_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_z_xyyz_zzz = cbuffer.data(gf_geom_01_off + 379 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyz_xxx, g_0_z_xyyz_xxy, g_0_z_xyyz_xxz, g_0_z_xyyz_xyy, g_0_z_xyyz_xyz, g_0_z_xyyz_xzz, g_0_z_xyyz_yyy, g_0_z_xyyz_yyz, g_0_z_xyyz_yzz, g_0_z_xyyz_zzz, g_0_z_yyz_xxx, g_0_z_yyz_xxxx, g_0_z_yyz_xxxy, g_0_z_yyz_xxxz, g_0_z_yyz_xxy, g_0_z_yyz_xxyy, g_0_z_yyz_xxyz, g_0_z_yyz_xxz, g_0_z_yyz_xxzz, g_0_z_yyz_xyy, g_0_z_yyz_xyyy, g_0_z_yyz_xyyz, g_0_z_yyz_xyz, g_0_z_yyz_xyzz, g_0_z_yyz_xzz, g_0_z_yyz_xzzz, g_0_z_yyz_yyy, g_0_z_yyz_yyz, g_0_z_yyz_yzz, g_0_z_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyz_xxx[k] = -g_0_z_yyz_xxx[k] * ab_x + g_0_z_yyz_xxxx[k];

                g_0_z_xyyz_xxy[k] = -g_0_z_yyz_xxy[k] * ab_x + g_0_z_yyz_xxxy[k];

                g_0_z_xyyz_xxz[k] = -g_0_z_yyz_xxz[k] * ab_x + g_0_z_yyz_xxxz[k];

                g_0_z_xyyz_xyy[k] = -g_0_z_yyz_xyy[k] * ab_x + g_0_z_yyz_xxyy[k];

                g_0_z_xyyz_xyz[k] = -g_0_z_yyz_xyz[k] * ab_x + g_0_z_yyz_xxyz[k];

                g_0_z_xyyz_xzz[k] = -g_0_z_yyz_xzz[k] * ab_x + g_0_z_yyz_xxzz[k];

                g_0_z_xyyz_yyy[k] = -g_0_z_yyz_yyy[k] * ab_x + g_0_z_yyz_xyyy[k];

                g_0_z_xyyz_yyz[k] = -g_0_z_yyz_yyz[k] * ab_x + g_0_z_yyz_xyyz[k];

                g_0_z_xyyz_yzz[k] = -g_0_z_yyz_yzz[k] * ab_x + g_0_z_yyz_xyzz[k];

                g_0_z_xyyz_zzz[k] = -g_0_z_yyz_zzz[k] * ab_x + g_0_z_yyz_xzzz[k];
            }

            /// Set up 380-390 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzz_xxx = cbuffer.data(gf_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_z_xyzz_xxy = cbuffer.data(gf_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_z_xyzz_xxz = cbuffer.data(gf_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_z_xyzz_xyy = cbuffer.data(gf_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_z_xyzz_xyz = cbuffer.data(gf_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_z_xyzz_xzz = cbuffer.data(gf_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_z_xyzz_yyy = cbuffer.data(gf_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_z_xyzz_yyz = cbuffer.data(gf_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_z_xyzz_yzz = cbuffer.data(gf_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_z_xyzz_zzz = cbuffer.data(gf_geom_01_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzz_xxx, g_0_z_xyzz_xxy, g_0_z_xyzz_xxz, g_0_z_xyzz_xyy, g_0_z_xyzz_xyz, g_0_z_xyzz_xzz, g_0_z_xyzz_yyy, g_0_z_xyzz_yyz, g_0_z_xyzz_yzz, g_0_z_xyzz_zzz, g_0_z_yzz_xxx, g_0_z_yzz_xxxx, g_0_z_yzz_xxxy, g_0_z_yzz_xxxz, g_0_z_yzz_xxy, g_0_z_yzz_xxyy, g_0_z_yzz_xxyz, g_0_z_yzz_xxz, g_0_z_yzz_xxzz, g_0_z_yzz_xyy, g_0_z_yzz_xyyy, g_0_z_yzz_xyyz, g_0_z_yzz_xyz, g_0_z_yzz_xyzz, g_0_z_yzz_xzz, g_0_z_yzz_xzzz, g_0_z_yzz_yyy, g_0_z_yzz_yyz, g_0_z_yzz_yzz, g_0_z_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzz_xxx[k] = -g_0_z_yzz_xxx[k] * ab_x + g_0_z_yzz_xxxx[k];

                g_0_z_xyzz_xxy[k] = -g_0_z_yzz_xxy[k] * ab_x + g_0_z_yzz_xxxy[k];

                g_0_z_xyzz_xxz[k] = -g_0_z_yzz_xxz[k] * ab_x + g_0_z_yzz_xxxz[k];

                g_0_z_xyzz_xyy[k] = -g_0_z_yzz_xyy[k] * ab_x + g_0_z_yzz_xxyy[k];

                g_0_z_xyzz_xyz[k] = -g_0_z_yzz_xyz[k] * ab_x + g_0_z_yzz_xxyz[k];

                g_0_z_xyzz_xzz[k] = -g_0_z_yzz_xzz[k] * ab_x + g_0_z_yzz_xxzz[k];

                g_0_z_xyzz_yyy[k] = -g_0_z_yzz_yyy[k] * ab_x + g_0_z_yzz_xyyy[k];

                g_0_z_xyzz_yyz[k] = -g_0_z_yzz_yyz[k] * ab_x + g_0_z_yzz_xyyz[k];

                g_0_z_xyzz_yzz[k] = -g_0_z_yzz_yzz[k] * ab_x + g_0_z_yzz_xyzz[k];

                g_0_z_xyzz_zzz[k] = -g_0_z_yzz_zzz[k] * ab_x + g_0_z_yzz_xzzz[k];
            }

            /// Set up 390-400 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzz_xxx = cbuffer.data(gf_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_z_xzzz_xxy = cbuffer.data(gf_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_z_xzzz_xxz = cbuffer.data(gf_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_z_xzzz_xyy = cbuffer.data(gf_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_z_xzzz_xyz = cbuffer.data(gf_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_z_xzzz_xzz = cbuffer.data(gf_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_z_xzzz_yyy = cbuffer.data(gf_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_z_xzzz_yyz = cbuffer.data(gf_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_z_xzzz_yzz = cbuffer.data(gf_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_z_xzzz_zzz = cbuffer.data(gf_geom_01_off + 399 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzz_xxx, g_0_z_xzzz_xxy, g_0_z_xzzz_xxz, g_0_z_xzzz_xyy, g_0_z_xzzz_xyz, g_0_z_xzzz_xzz, g_0_z_xzzz_yyy, g_0_z_xzzz_yyz, g_0_z_xzzz_yzz, g_0_z_xzzz_zzz, g_0_z_zzz_xxx, g_0_z_zzz_xxxx, g_0_z_zzz_xxxy, g_0_z_zzz_xxxz, g_0_z_zzz_xxy, g_0_z_zzz_xxyy, g_0_z_zzz_xxyz, g_0_z_zzz_xxz, g_0_z_zzz_xxzz, g_0_z_zzz_xyy, g_0_z_zzz_xyyy, g_0_z_zzz_xyyz, g_0_z_zzz_xyz, g_0_z_zzz_xyzz, g_0_z_zzz_xzz, g_0_z_zzz_xzzz, g_0_z_zzz_yyy, g_0_z_zzz_yyz, g_0_z_zzz_yzz, g_0_z_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzz_xxx[k] = -g_0_z_zzz_xxx[k] * ab_x + g_0_z_zzz_xxxx[k];

                g_0_z_xzzz_xxy[k] = -g_0_z_zzz_xxy[k] * ab_x + g_0_z_zzz_xxxy[k];

                g_0_z_xzzz_xxz[k] = -g_0_z_zzz_xxz[k] * ab_x + g_0_z_zzz_xxxz[k];

                g_0_z_xzzz_xyy[k] = -g_0_z_zzz_xyy[k] * ab_x + g_0_z_zzz_xxyy[k];

                g_0_z_xzzz_xyz[k] = -g_0_z_zzz_xyz[k] * ab_x + g_0_z_zzz_xxyz[k];

                g_0_z_xzzz_xzz[k] = -g_0_z_zzz_xzz[k] * ab_x + g_0_z_zzz_xxzz[k];

                g_0_z_xzzz_yyy[k] = -g_0_z_zzz_yyy[k] * ab_x + g_0_z_zzz_xyyy[k];

                g_0_z_xzzz_yyz[k] = -g_0_z_zzz_yyz[k] * ab_x + g_0_z_zzz_xyyz[k];

                g_0_z_xzzz_yzz[k] = -g_0_z_zzz_yzz[k] * ab_x + g_0_z_zzz_xyzz[k];

                g_0_z_xzzz_zzz[k] = -g_0_z_zzz_zzz[k] * ab_x + g_0_z_zzz_xzzz[k];
            }

            /// Set up 400-410 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyy_xxx = cbuffer.data(gf_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_z_yyyy_xxy = cbuffer.data(gf_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_z_yyyy_xxz = cbuffer.data(gf_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_z_yyyy_xyy = cbuffer.data(gf_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_z_yyyy_xyz = cbuffer.data(gf_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_z_yyyy_xzz = cbuffer.data(gf_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_z_yyyy_yyy = cbuffer.data(gf_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_z_yyyy_yyz = cbuffer.data(gf_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_z_yyyy_yzz = cbuffer.data(gf_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_z_yyyy_zzz = cbuffer.data(gf_geom_01_off + 409 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyy_xxx, g_0_z_yyy_xxxy, g_0_z_yyy_xxy, g_0_z_yyy_xxyy, g_0_z_yyy_xxyz, g_0_z_yyy_xxz, g_0_z_yyy_xyy, g_0_z_yyy_xyyy, g_0_z_yyy_xyyz, g_0_z_yyy_xyz, g_0_z_yyy_xyzz, g_0_z_yyy_xzz, g_0_z_yyy_yyy, g_0_z_yyy_yyyy, g_0_z_yyy_yyyz, g_0_z_yyy_yyz, g_0_z_yyy_yyzz, g_0_z_yyy_yzz, g_0_z_yyy_yzzz, g_0_z_yyy_zzz, g_0_z_yyyy_xxx, g_0_z_yyyy_xxy, g_0_z_yyyy_xxz, g_0_z_yyyy_xyy, g_0_z_yyyy_xyz, g_0_z_yyyy_xzz, g_0_z_yyyy_yyy, g_0_z_yyyy_yyz, g_0_z_yyyy_yzz, g_0_z_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyy_xxx[k] = -g_0_z_yyy_xxx[k] * ab_y + g_0_z_yyy_xxxy[k];

                g_0_z_yyyy_xxy[k] = -g_0_z_yyy_xxy[k] * ab_y + g_0_z_yyy_xxyy[k];

                g_0_z_yyyy_xxz[k] = -g_0_z_yyy_xxz[k] * ab_y + g_0_z_yyy_xxyz[k];

                g_0_z_yyyy_xyy[k] = -g_0_z_yyy_xyy[k] * ab_y + g_0_z_yyy_xyyy[k];

                g_0_z_yyyy_xyz[k] = -g_0_z_yyy_xyz[k] * ab_y + g_0_z_yyy_xyyz[k];

                g_0_z_yyyy_xzz[k] = -g_0_z_yyy_xzz[k] * ab_y + g_0_z_yyy_xyzz[k];

                g_0_z_yyyy_yyy[k] = -g_0_z_yyy_yyy[k] * ab_y + g_0_z_yyy_yyyy[k];

                g_0_z_yyyy_yyz[k] = -g_0_z_yyy_yyz[k] * ab_y + g_0_z_yyy_yyyz[k];

                g_0_z_yyyy_yzz[k] = -g_0_z_yyy_yzz[k] * ab_y + g_0_z_yyy_yyzz[k];

                g_0_z_yyyy_zzz[k] = -g_0_z_yyy_zzz[k] * ab_y + g_0_z_yyy_yzzz[k];
            }

            /// Set up 410-420 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyz_xxx = cbuffer.data(gf_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_z_yyyz_xxy = cbuffer.data(gf_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_z_yyyz_xxz = cbuffer.data(gf_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_z_yyyz_xyy = cbuffer.data(gf_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_z_yyyz_xyz = cbuffer.data(gf_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_z_yyyz_xzz = cbuffer.data(gf_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_z_yyyz_yyy = cbuffer.data(gf_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_z_yyyz_yyz = cbuffer.data(gf_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_z_yyyz_yzz = cbuffer.data(gf_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_z_yyyz_zzz = cbuffer.data(gf_geom_01_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyz_xxx, g_0_z_yyyz_xxy, g_0_z_yyyz_xxz, g_0_z_yyyz_xyy, g_0_z_yyyz_xyz, g_0_z_yyyz_xzz, g_0_z_yyyz_yyy, g_0_z_yyyz_yyz, g_0_z_yyyz_yzz, g_0_z_yyyz_zzz, g_0_z_yyz_xxx, g_0_z_yyz_xxxy, g_0_z_yyz_xxy, g_0_z_yyz_xxyy, g_0_z_yyz_xxyz, g_0_z_yyz_xxz, g_0_z_yyz_xyy, g_0_z_yyz_xyyy, g_0_z_yyz_xyyz, g_0_z_yyz_xyz, g_0_z_yyz_xyzz, g_0_z_yyz_xzz, g_0_z_yyz_yyy, g_0_z_yyz_yyyy, g_0_z_yyz_yyyz, g_0_z_yyz_yyz, g_0_z_yyz_yyzz, g_0_z_yyz_yzz, g_0_z_yyz_yzzz, g_0_z_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyz_xxx[k] = -g_0_z_yyz_xxx[k] * ab_y + g_0_z_yyz_xxxy[k];

                g_0_z_yyyz_xxy[k] = -g_0_z_yyz_xxy[k] * ab_y + g_0_z_yyz_xxyy[k];

                g_0_z_yyyz_xxz[k] = -g_0_z_yyz_xxz[k] * ab_y + g_0_z_yyz_xxyz[k];

                g_0_z_yyyz_xyy[k] = -g_0_z_yyz_xyy[k] * ab_y + g_0_z_yyz_xyyy[k];

                g_0_z_yyyz_xyz[k] = -g_0_z_yyz_xyz[k] * ab_y + g_0_z_yyz_xyyz[k];

                g_0_z_yyyz_xzz[k] = -g_0_z_yyz_xzz[k] * ab_y + g_0_z_yyz_xyzz[k];

                g_0_z_yyyz_yyy[k] = -g_0_z_yyz_yyy[k] * ab_y + g_0_z_yyz_yyyy[k];

                g_0_z_yyyz_yyz[k] = -g_0_z_yyz_yyz[k] * ab_y + g_0_z_yyz_yyyz[k];

                g_0_z_yyyz_yzz[k] = -g_0_z_yyz_yzz[k] * ab_y + g_0_z_yyz_yyzz[k];

                g_0_z_yyyz_zzz[k] = -g_0_z_yyz_zzz[k] * ab_y + g_0_z_yyz_yzzz[k];
            }

            /// Set up 420-430 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzz_xxx = cbuffer.data(gf_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_yyzz_xxy = cbuffer.data(gf_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_yyzz_xxz = cbuffer.data(gf_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_yyzz_xyy = cbuffer.data(gf_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_yyzz_xyz = cbuffer.data(gf_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_yyzz_xzz = cbuffer.data(gf_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_z_yyzz_yyy = cbuffer.data(gf_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_yyzz_yyz = cbuffer.data(gf_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_yyzz_yzz = cbuffer.data(gf_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_yyzz_zzz = cbuffer.data(gf_geom_01_off + 429 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzz_xxx, g_0_z_yyzz_xxy, g_0_z_yyzz_xxz, g_0_z_yyzz_xyy, g_0_z_yyzz_xyz, g_0_z_yyzz_xzz, g_0_z_yyzz_yyy, g_0_z_yyzz_yyz, g_0_z_yyzz_yzz, g_0_z_yyzz_zzz, g_0_z_yzz_xxx, g_0_z_yzz_xxxy, g_0_z_yzz_xxy, g_0_z_yzz_xxyy, g_0_z_yzz_xxyz, g_0_z_yzz_xxz, g_0_z_yzz_xyy, g_0_z_yzz_xyyy, g_0_z_yzz_xyyz, g_0_z_yzz_xyz, g_0_z_yzz_xyzz, g_0_z_yzz_xzz, g_0_z_yzz_yyy, g_0_z_yzz_yyyy, g_0_z_yzz_yyyz, g_0_z_yzz_yyz, g_0_z_yzz_yyzz, g_0_z_yzz_yzz, g_0_z_yzz_yzzz, g_0_z_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzz_xxx[k] = -g_0_z_yzz_xxx[k] * ab_y + g_0_z_yzz_xxxy[k];

                g_0_z_yyzz_xxy[k] = -g_0_z_yzz_xxy[k] * ab_y + g_0_z_yzz_xxyy[k];

                g_0_z_yyzz_xxz[k] = -g_0_z_yzz_xxz[k] * ab_y + g_0_z_yzz_xxyz[k];

                g_0_z_yyzz_xyy[k] = -g_0_z_yzz_xyy[k] * ab_y + g_0_z_yzz_xyyy[k];

                g_0_z_yyzz_xyz[k] = -g_0_z_yzz_xyz[k] * ab_y + g_0_z_yzz_xyyz[k];

                g_0_z_yyzz_xzz[k] = -g_0_z_yzz_xzz[k] * ab_y + g_0_z_yzz_xyzz[k];

                g_0_z_yyzz_yyy[k] = -g_0_z_yzz_yyy[k] * ab_y + g_0_z_yzz_yyyy[k];

                g_0_z_yyzz_yyz[k] = -g_0_z_yzz_yyz[k] * ab_y + g_0_z_yzz_yyyz[k];

                g_0_z_yyzz_yzz[k] = -g_0_z_yzz_yzz[k] * ab_y + g_0_z_yzz_yyzz[k];

                g_0_z_yyzz_zzz[k] = -g_0_z_yzz_zzz[k] * ab_y + g_0_z_yzz_yzzz[k];
            }

            /// Set up 430-440 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzz_xxx = cbuffer.data(gf_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_yzzz_xxy = cbuffer.data(gf_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_yzzz_xxz = cbuffer.data(gf_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_yzzz_xyy = cbuffer.data(gf_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_yzzz_xyz = cbuffer.data(gf_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_yzzz_xzz = cbuffer.data(gf_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_yzzz_yyy = cbuffer.data(gf_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_yzzz_yyz = cbuffer.data(gf_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_yzzz_yzz = cbuffer.data(gf_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_yzzz_zzz = cbuffer.data(gf_geom_01_off + 439 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzz_xxx, g_0_z_yzzz_xxy, g_0_z_yzzz_xxz, g_0_z_yzzz_xyy, g_0_z_yzzz_xyz, g_0_z_yzzz_xzz, g_0_z_yzzz_yyy, g_0_z_yzzz_yyz, g_0_z_yzzz_yzz, g_0_z_yzzz_zzz, g_0_z_zzz_xxx, g_0_z_zzz_xxxy, g_0_z_zzz_xxy, g_0_z_zzz_xxyy, g_0_z_zzz_xxyz, g_0_z_zzz_xxz, g_0_z_zzz_xyy, g_0_z_zzz_xyyy, g_0_z_zzz_xyyz, g_0_z_zzz_xyz, g_0_z_zzz_xyzz, g_0_z_zzz_xzz, g_0_z_zzz_yyy, g_0_z_zzz_yyyy, g_0_z_zzz_yyyz, g_0_z_zzz_yyz, g_0_z_zzz_yyzz, g_0_z_zzz_yzz, g_0_z_zzz_yzzz, g_0_z_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzz_xxx[k] = -g_0_z_zzz_xxx[k] * ab_y + g_0_z_zzz_xxxy[k];

                g_0_z_yzzz_xxy[k] = -g_0_z_zzz_xxy[k] * ab_y + g_0_z_zzz_xxyy[k];

                g_0_z_yzzz_xxz[k] = -g_0_z_zzz_xxz[k] * ab_y + g_0_z_zzz_xxyz[k];

                g_0_z_yzzz_xyy[k] = -g_0_z_zzz_xyy[k] * ab_y + g_0_z_zzz_xyyy[k];

                g_0_z_yzzz_xyz[k] = -g_0_z_zzz_xyz[k] * ab_y + g_0_z_zzz_xyyz[k];

                g_0_z_yzzz_xzz[k] = -g_0_z_zzz_xzz[k] * ab_y + g_0_z_zzz_xyzz[k];

                g_0_z_yzzz_yyy[k] = -g_0_z_zzz_yyy[k] * ab_y + g_0_z_zzz_yyyy[k];

                g_0_z_yzzz_yyz[k] = -g_0_z_zzz_yyz[k] * ab_y + g_0_z_zzz_yyyz[k];

                g_0_z_yzzz_yzz[k] = -g_0_z_zzz_yzz[k] * ab_y + g_0_z_zzz_yyzz[k];

                g_0_z_yzzz_zzz[k] = -g_0_z_zzz_zzz[k] * ab_y + g_0_z_zzz_yzzz[k];
            }

            /// Set up 440-450 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzz_xxx = cbuffer.data(gf_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_zzzz_xxy = cbuffer.data(gf_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_zzzz_xxz = cbuffer.data(gf_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_zzzz_xyy = cbuffer.data(gf_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_zzzz_xyz = cbuffer.data(gf_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_zzzz_xzz = cbuffer.data(gf_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_zzzz_yyy = cbuffer.data(gf_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_zzzz_yyz = cbuffer.data(gf_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_zzzz_yzz = cbuffer.data(gf_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_zzzz_zzz = cbuffer.data(gf_geom_01_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_xxx, g_0_z_zzz_xxxz, g_0_z_zzz_xxy, g_0_z_zzz_xxyz, g_0_z_zzz_xxz, g_0_z_zzz_xxzz, g_0_z_zzz_xyy, g_0_z_zzz_xyyz, g_0_z_zzz_xyz, g_0_z_zzz_xyzz, g_0_z_zzz_xzz, g_0_z_zzz_xzzz, g_0_z_zzz_yyy, g_0_z_zzz_yyyz, g_0_z_zzz_yyz, g_0_z_zzz_yyzz, g_0_z_zzz_yzz, g_0_z_zzz_yzzz, g_0_z_zzz_zzz, g_0_z_zzz_zzzz, g_0_z_zzzz_xxx, g_0_z_zzzz_xxy, g_0_z_zzzz_xxz, g_0_z_zzzz_xyy, g_0_z_zzzz_xyz, g_0_z_zzzz_xzz, g_0_z_zzzz_yyy, g_0_z_zzzz_yyz, g_0_z_zzzz_yzz, g_0_z_zzzz_zzz, g_zzz_xxx, g_zzz_xxy, g_zzz_xxz, g_zzz_xyy, g_zzz_xyz, g_zzz_xzz, g_zzz_yyy, g_zzz_yyz, g_zzz_yzz, g_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzz_xxx[k] = g_zzz_xxx[k] - g_0_z_zzz_xxx[k] * ab_z + g_0_z_zzz_xxxz[k];

                g_0_z_zzzz_xxy[k] = g_zzz_xxy[k] - g_0_z_zzz_xxy[k] * ab_z + g_0_z_zzz_xxyz[k];

                g_0_z_zzzz_xxz[k] = g_zzz_xxz[k] - g_0_z_zzz_xxz[k] * ab_z + g_0_z_zzz_xxzz[k];

                g_0_z_zzzz_xyy[k] = g_zzz_xyy[k] - g_0_z_zzz_xyy[k] * ab_z + g_0_z_zzz_xyyz[k];

                g_0_z_zzzz_xyz[k] = g_zzz_xyz[k] - g_0_z_zzz_xyz[k] * ab_z + g_0_z_zzz_xyzz[k];

                g_0_z_zzzz_xzz[k] = g_zzz_xzz[k] - g_0_z_zzz_xzz[k] * ab_z + g_0_z_zzz_xzzz[k];

                g_0_z_zzzz_yyy[k] = g_zzz_yyy[k] - g_0_z_zzz_yyy[k] * ab_z + g_0_z_zzz_yyyz[k];

                g_0_z_zzzz_yyz[k] = g_zzz_yyz[k] - g_0_z_zzz_yyz[k] * ab_z + g_0_z_zzz_yyzz[k];

                g_0_z_zzzz_yzz[k] = g_zzz_yzz[k] - g_0_z_zzz_yzz[k] * ab_z + g_0_z_zzz_yzzz[k];

                g_0_z_zzzz_zzz[k] = g_zzz_zzz[k] - g_0_z_zzz_zzz[k] * ab_z + g_0_z_zzz_zzzz[k];
            }
        }
    }
}

} // erirec namespace

