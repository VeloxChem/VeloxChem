#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_pdxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_pdxx,
                                              const size_t idx_geom_0010_sdxx,
                                              const size_t idx_geom_1010_sdxx,
                                              const size_t idx_geom_1010_sfxx,
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
            /// Set up components of auxilary buffer : SDSS

            const auto sd_geom_0010_off = idx_geom_0010_sdxx + i * dcomps + j;

            auto g_0_0_x_0_0_xx = cbuffer.data(sd_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_0_xy = cbuffer.data(sd_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_0_xz = cbuffer.data(sd_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_0_yy = cbuffer.data(sd_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_0_yz = cbuffer.data(sd_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_0_zz = cbuffer.data(sd_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_y_0_0_xx = cbuffer.data(sd_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_y_0_0_xy = cbuffer.data(sd_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_y_0_0_xz = cbuffer.data(sd_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_y_0_0_yy = cbuffer.data(sd_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_y_0_0_yz = cbuffer.data(sd_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_y_0_0_zz = cbuffer.data(sd_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_z_0_0_xx = cbuffer.data(sd_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_z_0_0_xy = cbuffer.data(sd_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_z_0_0_xz = cbuffer.data(sd_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_z_0_0_yy = cbuffer.data(sd_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_z_0_0_yz = cbuffer.data(sd_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_z_0_0_zz = cbuffer.data(sd_geom_0010_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SDSS

            const auto sd_geom_1010_off = idx_geom_1010_sdxx + i * dcomps + j;

            auto g_x_0_x_0_0_xx = cbuffer.data(sd_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_0_xy = cbuffer.data(sd_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_0_xz = cbuffer.data(sd_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_0_yy = cbuffer.data(sd_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_0_yz = cbuffer.data(sd_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_0_zz = cbuffer.data(sd_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_y_0_0_xx = cbuffer.data(sd_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_y_0_0_xy = cbuffer.data(sd_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_y_0_0_xz = cbuffer.data(sd_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_y_0_0_yy = cbuffer.data(sd_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_y_0_0_yz = cbuffer.data(sd_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_y_0_0_zz = cbuffer.data(sd_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_z_0_0_xx = cbuffer.data(sd_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_z_0_0_xy = cbuffer.data(sd_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_z_0_0_xz = cbuffer.data(sd_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_z_0_0_yy = cbuffer.data(sd_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_z_0_0_yz = cbuffer.data(sd_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_z_0_0_zz = cbuffer.data(sd_geom_1010_off + 17 * ccomps * dcomps);

            auto g_y_0_x_0_0_xx = cbuffer.data(sd_geom_1010_off + 18 * ccomps * dcomps);

            auto g_y_0_x_0_0_xy = cbuffer.data(sd_geom_1010_off + 19 * ccomps * dcomps);

            auto g_y_0_x_0_0_xz = cbuffer.data(sd_geom_1010_off + 20 * ccomps * dcomps);

            auto g_y_0_x_0_0_yy = cbuffer.data(sd_geom_1010_off + 21 * ccomps * dcomps);

            auto g_y_0_x_0_0_yz = cbuffer.data(sd_geom_1010_off + 22 * ccomps * dcomps);

            auto g_y_0_x_0_0_zz = cbuffer.data(sd_geom_1010_off + 23 * ccomps * dcomps);

            auto g_y_0_y_0_0_xx = cbuffer.data(sd_geom_1010_off + 24 * ccomps * dcomps);

            auto g_y_0_y_0_0_xy = cbuffer.data(sd_geom_1010_off + 25 * ccomps * dcomps);

            auto g_y_0_y_0_0_xz = cbuffer.data(sd_geom_1010_off + 26 * ccomps * dcomps);

            auto g_y_0_y_0_0_yy = cbuffer.data(sd_geom_1010_off + 27 * ccomps * dcomps);

            auto g_y_0_y_0_0_yz = cbuffer.data(sd_geom_1010_off + 28 * ccomps * dcomps);

            auto g_y_0_y_0_0_zz = cbuffer.data(sd_geom_1010_off + 29 * ccomps * dcomps);

            auto g_y_0_z_0_0_xx = cbuffer.data(sd_geom_1010_off + 30 * ccomps * dcomps);

            auto g_y_0_z_0_0_xy = cbuffer.data(sd_geom_1010_off + 31 * ccomps * dcomps);

            auto g_y_0_z_0_0_xz = cbuffer.data(sd_geom_1010_off + 32 * ccomps * dcomps);

            auto g_y_0_z_0_0_yy = cbuffer.data(sd_geom_1010_off + 33 * ccomps * dcomps);

            auto g_y_0_z_0_0_yz = cbuffer.data(sd_geom_1010_off + 34 * ccomps * dcomps);

            auto g_y_0_z_0_0_zz = cbuffer.data(sd_geom_1010_off + 35 * ccomps * dcomps);

            auto g_z_0_x_0_0_xx = cbuffer.data(sd_geom_1010_off + 36 * ccomps * dcomps);

            auto g_z_0_x_0_0_xy = cbuffer.data(sd_geom_1010_off + 37 * ccomps * dcomps);

            auto g_z_0_x_0_0_xz = cbuffer.data(sd_geom_1010_off + 38 * ccomps * dcomps);

            auto g_z_0_x_0_0_yy = cbuffer.data(sd_geom_1010_off + 39 * ccomps * dcomps);

            auto g_z_0_x_0_0_yz = cbuffer.data(sd_geom_1010_off + 40 * ccomps * dcomps);

            auto g_z_0_x_0_0_zz = cbuffer.data(sd_geom_1010_off + 41 * ccomps * dcomps);

            auto g_z_0_y_0_0_xx = cbuffer.data(sd_geom_1010_off + 42 * ccomps * dcomps);

            auto g_z_0_y_0_0_xy = cbuffer.data(sd_geom_1010_off + 43 * ccomps * dcomps);

            auto g_z_0_y_0_0_xz = cbuffer.data(sd_geom_1010_off + 44 * ccomps * dcomps);

            auto g_z_0_y_0_0_yy = cbuffer.data(sd_geom_1010_off + 45 * ccomps * dcomps);

            auto g_z_0_y_0_0_yz = cbuffer.data(sd_geom_1010_off + 46 * ccomps * dcomps);

            auto g_z_0_y_0_0_zz = cbuffer.data(sd_geom_1010_off + 47 * ccomps * dcomps);

            auto g_z_0_z_0_0_xx = cbuffer.data(sd_geom_1010_off + 48 * ccomps * dcomps);

            auto g_z_0_z_0_0_xy = cbuffer.data(sd_geom_1010_off + 49 * ccomps * dcomps);

            auto g_z_0_z_0_0_xz = cbuffer.data(sd_geom_1010_off + 50 * ccomps * dcomps);

            auto g_z_0_z_0_0_yy = cbuffer.data(sd_geom_1010_off + 51 * ccomps * dcomps);

            auto g_z_0_z_0_0_yz = cbuffer.data(sd_geom_1010_off + 52 * ccomps * dcomps);

            auto g_z_0_z_0_0_zz = cbuffer.data(sd_geom_1010_off + 53 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SFSS

            const auto sf_geom_1010_off = idx_geom_1010_sfxx + i * dcomps + j;

            auto g_x_0_x_0_0_xxx = cbuffer.data(sf_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxy = cbuffer.data(sf_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxz = cbuffer.data(sf_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyy = cbuffer.data(sf_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyz = cbuffer.data(sf_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_0_xzz = cbuffer.data(sf_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyy = cbuffer.data(sf_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyz = cbuffer.data(sf_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_0_yzz = cbuffer.data(sf_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_0_zzz = cbuffer.data(sf_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxx = cbuffer.data(sf_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxy = cbuffer.data(sf_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxz = cbuffer.data(sf_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyy = cbuffer.data(sf_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyz = cbuffer.data(sf_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_y_0_0_xzz = cbuffer.data(sf_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyy = cbuffer.data(sf_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyz = cbuffer.data(sf_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_y_0_0_yzz = cbuffer.data(sf_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_y_0_0_zzz = cbuffer.data(sf_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxx = cbuffer.data(sf_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxy = cbuffer.data(sf_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxz = cbuffer.data(sf_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyy = cbuffer.data(sf_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyz = cbuffer.data(sf_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_z_0_0_xzz = cbuffer.data(sf_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyy = cbuffer.data(sf_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyz = cbuffer.data(sf_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_z_0_0_yzz = cbuffer.data(sf_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_z_0_0_zzz = cbuffer.data(sf_geom_1010_off + 29 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxx = cbuffer.data(sf_geom_1010_off + 30 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxy = cbuffer.data(sf_geom_1010_off + 31 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxz = cbuffer.data(sf_geom_1010_off + 32 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyy = cbuffer.data(sf_geom_1010_off + 33 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyz = cbuffer.data(sf_geom_1010_off + 34 * ccomps * dcomps);

            auto g_y_0_x_0_0_xzz = cbuffer.data(sf_geom_1010_off + 35 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyy = cbuffer.data(sf_geom_1010_off + 36 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyz = cbuffer.data(sf_geom_1010_off + 37 * ccomps * dcomps);

            auto g_y_0_x_0_0_yzz = cbuffer.data(sf_geom_1010_off + 38 * ccomps * dcomps);

            auto g_y_0_x_0_0_zzz = cbuffer.data(sf_geom_1010_off + 39 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxx = cbuffer.data(sf_geom_1010_off + 40 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxy = cbuffer.data(sf_geom_1010_off + 41 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxz = cbuffer.data(sf_geom_1010_off + 42 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyy = cbuffer.data(sf_geom_1010_off + 43 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyz = cbuffer.data(sf_geom_1010_off + 44 * ccomps * dcomps);

            auto g_y_0_y_0_0_xzz = cbuffer.data(sf_geom_1010_off + 45 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyy = cbuffer.data(sf_geom_1010_off + 46 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyz = cbuffer.data(sf_geom_1010_off + 47 * ccomps * dcomps);

            auto g_y_0_y_0_0_yzz = cbuffer.data(sf_geom_1010_off + 48 * ccomps * dcomps);

            auto g_y_0_y_0_0_zzz = cbuffer.data(sf_geom_1010_off + 49 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxx = cbuffer.data(sf_geom_1010_off + 50 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxy = cbuffer.data(sf_geom_1010_off + 51 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxz = cbuffer.data(sf_geom_1010_off + 52 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyy = cbuffer.data(sf_geom_1010_off + 53 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyz = cbuffer.data(sf_geom_1010_off + 54 * ccomps * dcomps);

            auto g_y_0_z_0_0_xzz = cbuffer.data(sf_geom_1010_off + 55 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyy = cbuffer.data(sf_geom_1010_off + 56 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyz = cbuffer.data(sf_geom_1010_off + 57 * ccomps * dcomps);

            auto g_y_0_z_0_0_yzz = cbuffer.data(sf_geom_1010_off + 58 * ccomps * dcomps);

            auto g_y_0_z_0_0_zzz = cbuffer.data(sf_geom_1010_off + 59 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxx = cbuffer.data(sf_geom_1010_off + 60 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxy = cbuffer.data(sf_geom_1010_off + 61 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxz = cbuffer.data(sf_geom_1010_off + 62 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyy = cbuffer.data(sf_geom_1010_off + 63 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyz = cbuffer.data(sf_geom_1010_off + 64 * ccomps * dcomps);

            auto g_z_0_x_0_0_xzz = cbuffer.data(sf_geom_1010_off + 65 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyy = cbuffer.data(sf_geom_1010_off + 66 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyz = cbuffer.data(sf_geom_1010_off + 67 * ccomps * dcomps);

            auto g_z_0_x_0_0_yzz = cbuffer.data(sf_geom_1010_off + 68 * ccomps * dcomps);

            auto g_z_0_x_0_0_zzz = cbuffer.data(sf_geom_1010_off + 69 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxx = cbuffer.data(sf_geom_1010_off + 70 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxy = cbuffer.data(sf_geom_1010_off + 71 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxz = cbuffer.data(sf_geom_1010_off + 72 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyy = cbuffer.data(sf_geom_1010_off + 73 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyz = cbuffer.data(sf_geom_1010_off + 74 * ccomps * dcomps);

            auto g_z_0_y_0_0_xzz = cbuffer.data(sf_geom_1010_off + 75 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyy = cbuffer.data(sf_geom_1010_off + 76 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyz = cbuffer.data(sf_geom_1010_off + 77 * ccomps * dcomps);

            auto g_z_0_y_0_0_yzz = cbuffer.data(sf_geom_1010_off + 78 * ccomps * dcomps);

            auto g_z_0_y_0_0_zzz = cbuffer.data(sf_geom_1010_off + 79 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxx = cbuffer.data(sf_geom_1010_off + 80 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxy = cbuffer.data(sf_geom_1010_off + 81 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxz = cbuffer.data(sf_geom_1010_off + 82 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyy = cbuffer.data(sf_geom_1010_off + 83 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyz = cbuffer.data(sf_geom_1010_off + 84 * ccomps * dcomps);

            auto g_z_0_z_0_0_xzz = cbuffer.data(sf_geom_1010_off + 85 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyy = cbuffer.data(sf_geom_1010_off + 86 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyz = cbuffer.data(sf_geom_1010_off + 87 * ccomps * dcomps);

            auto g_z_0_z_0_0_yzz = cbuffer.data(sf_geom_1010_off + 88 * ccomps * dcomps);

            auto g_z_0_z_0_0_zzz = cbuffer.data(sf_geom_1010_off + 89 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pdxx

            const auto pd_geom_1010_off = idx_geom_1010_pdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_x_xx = cbuffer.data(pd_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_x_xy = cbuffer.data(pd_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_x_xz = cbuffer.data(pd_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_x_yy = cbuffer.data(pd_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_x_yz = cbuffer.data(pd_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_x_zz = cbuffer.data(pd_geom_1010_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_xx, g_0_0_x_0_0_xy, g_0_0_x_0_0_xz, g_0_0_x_0_0_yy, g_0_0_x_0_0_yz, g_0_0_x_0_0_zz, g_x_0_x_0_0_xx, g_x_0_x_0_0_xxx, g_x_0_x_0_0_xxy, g_x_0_x_0_0_xxz, g_x_0_x_0_0_xy, g_x_0_x_0_0_xyy, g_x_0_x_0_0_xyz, g_x_0_x_0_0_xz, g_x_0_x_0_0_xzz, g_x_0_x_0_0_yy, g_x_0_x_0_0_yz, g_x_0_x_0_0_zz, g_x_0_x_0_x_xx, g_x_0_x_0_x_xy, g_x_0_x_0_x_xz, g_x_0_x_0_x_yy, g_x_0_x_0_x_yz, g_x_0_x_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_x_xx[k] = -g_0_0_x_0_0_xx[k] - g_x_0_x_0_0_xx[k] * ab_x + g_x_0_x_0_0_xxx[k];

                g_x_0_x_0_x_xy[k] = -g_0_0_x_0_0_xy[k] - g_x_0_x_0_0_xy[k] * ab_x + g_x_0_x_0_0_xxy[k];

                g_x_0_x_0_x_xz[k] = -g_0_0_x_0_0_xz[k] - g_x_0_x_0_0_xz[k] * ab_x + g_x_0_x_0_0_xxz[k];

                g_x_0_x_0_x_yy[k] = -g_0_0_x_0_0_yy[k] - g_x_0_x_0_0_yy[k] * ab_x + g_x_0_x_0_0_xyy[k];

                g_x_0_x_0_x_yz[k] = -g_0_0_x_0_0_yz[k] - g_x_0_x_0_0_yz[k] * ab_x + g_x_0_x_0_0_xyz[k];

                g_x_0_x_0_x_zz[k] = -g_0_0_x_0_0_zz[k] - g_x_0_x_0_0_zz[k] * ab_x + g_x_0_x_0_0_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_y_xx = cbuffer.data(pd_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_y_xy = cbuffer.data(pd_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_y_xz = cbuffer.data(pd_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_y_yy = cbuffer.data(pd_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_y_yz = cbuffer.data(pd_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_y_zz = cbuffer.data(pd_geom_1010_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_0_xx, g_x_0_x_0_0_xxy, g_x_0_x_0_0_xy, g_x_0_x_0_0_xyy, g_x_0_x_0_0_xyz, g_x_0_x_0_0_xz, g_x_0_x_0_0_yy, g_x_0_x_0_0_yyy, g_x_0_x_0_0_yyz, g_x_0_x_0_0_yz, g_x_0_x_0_0_yzz, g_x_0_x_0_0_zz, g_x_0_x_0_y_xx, g_x_0_x_0_y_xy, g_x_0_x_0_y_xz, g_x_0_x_0_y_yy, g_x_0_x_0_y_yz, g_x_0_x_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_y_xx[k] = -g_x_0_x_0_0_xx[k] * ab_y + g_x_0_x_0_0_xxy[k];

                g_x_0_x_0_y_xy[k] = -g_x_0_x_0_0_xy[k] * ab_y + g_x_0_x_0_0_xyy[k];

                g_x_0_x_0_y_xz[k] = -g_x_0_x_0_0_xz[k] * ab_y + g_x_0_x_0_0_xyz[k];

                g_x_0_x_0_y_yy[k] = -g_x_0_x_0_0_yy[k] * ab_y + g_x_0_x_0_0_yyy[k];

                g_x_0_x_0_y_yz[k] = -g_x_0_x_0_0_yz[k] * ab_y + g_x_0_x_0_0_yyz[k];

                g_x_0_x_0_y_zz[k] = -g_x_0_x_0_0_zz[k] * ab_y + g_x_0_x_0_0_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_z_xx = cbuffer.data(pd_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_z_xy = cbuffer.data(pd_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_z_xz = cbuffer.data(pd_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_z_yy = cbuffer.data(pd_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_z_yz = cbuffer.data(pd_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_z_zz = cbuffer.data(pd_geom_1010_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_0_xx, g_x_0_x_0_0_xxz, g_x_0_x_0_0_xy, g_x_0_x_0_0_xyz, g_x_0_x_0_0_xz, g_x_0_x_0_0_xzz, g_x_0_x_0_0_yy, g_x_0_x_0_0_yyz, g_x_0_x_0_0_yz, g_x_0_x_0_0_yzz, g_x_0_x_0_0_zz, g_x_0_x_0_0_zzz, g_x_0_x_0_z_xx, g_x_0_x_0_z_xy, g_x_0_x_0_z_xz, g_x_0_x_0_z_yy, g_x_0_x_0_z_yz, g_x_0_x_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_z_xx[k] = -g_x_0_x_0_0_xx[k] * ab_z + g_x_0_x_0_0_xxz[k];

                g_x_0_x_0_z_xy[k] = -g_x_0_x_0_0_xy[k] * ab_z + g_x_0_x_0_0_xyz[k];

                g_x_0_x_0_z_xz[k] = -g_x_0_x_0_0_xz[k] * ab_z + g_x_0_x_0_0_xzz[k];

                g_x_0_x_0_z_yy[k] = -g_x_0_x_0_0_yy[k] * ab_z + g_x_0_x_0_0_yyz[k];

                g_x_0_x_0_z_yz[k] = -g_x_0_x_0_0_yz[k] * ab_z + g_x_0_x_0_0_yzz[k];

                g_x_0_x_0_z_zz[k] = -g_x_0_x_0_0_zz[k] * ab_z + g_x_0_x_0_0_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_x_xx = cbuffer.data(pd_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_y_0_x_xy = cbuffer.data(pd_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_y_0_x_xz = cbuffer.data(pd_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_y_0_x_yy = cbuffer.data(pd_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_y_0_x_yz = cbuffer.data(pd_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_y_0_x_zz = cbuffer.data(pd_geom_1010_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_xx, g_0_0_y_0_0_xy, g_0_0_y_0_0_xz, g_0_0_y_0_0_yy, g_0_0_y_0_0_yz, g_0_0_y_0_0_zz, g_x_0_y_0_0_xx, g_x_0_y_0_0_xxx, g_x_0_y_0_0_xxy, g_x_0_y_0_0_xxz, g_x_0_y_0_0_xy, g_x_0_y_0_0_xyy, g_x_0_y_0_0_xyz, g_x_0_y_0_0_xz, g_x_0_y_0_0_xzz, g_x_0_y_0_0_yy, g_x_0_y_0_0_yz, g_x_0_y_0_0_zz, g_x_0_y_0_x_xx, g_x_0_y_0_x_xy, g_x_0_y_0_x_xz, g_x_0_y_0_x_yy, g_x_0_y_0_x_yz, g_x_0_y_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_x_xx[k] = -g_0_0_y_0_0_xx[k] - g_x_0_y_0_0_xx[k] * ab_x + g_x_0_y_0_0_xxx[k];

                g_x_0_y_0_x_xy[k] = -g_0_0_y_0_0_xy[k] - g_x_0_y_0_0_xy[k] * ab_x + g_x_0_y_0_0_xxy[k];

                g_x_0_y_0_x_xz[k] = -g_0_0_y_0_0_xz[k] - g_x_0_y_0_0_xz[k] * ab_x + g_x_0_y_0_0_xxz[k];

                g_x_0_y_0_x_yy[k] = -g_0_0_y_0_0_yy[k] - g_x_0_y_0_0_yy[k] * ab_x + g_x_0_y_0_0_xyy[k];

                g_x_0_y_0_x_yz[k] = -g_0_0_y_0_0_yz[k] - g_x_0_y_0_0_yz[k] * ab_x + g_x_0_y_0_0_xyz[k];

                g_x_0_y_0_x_zz[k] = -g_0_0_y_0_0_zz[k] - g_x_0_y_0_0_zz[k] * ab_x + g_x_0_y_0_0_xzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_y_xx = cbuffer.data(pd_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_y_0_y_xy = cbuffer.data(pd_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_y_0_y_xz = cbuffer.data(pd_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_y_0_y_yy = cbuffer.data(pd_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_y_0_y_yz = cbuffer.data(pd_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_y_0_y_zz = cbuffer.data(pd_geom_1010_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_0_xx, g_x_0_y_0_0_xxy, g_x_0_y_0_0_xy, g_x_0_y_0_0_xyy, g_x_0_y_0_0_xyz, g_x_0_y_0_0_xz, g_x_0_y_0_0_yy, g_x_0_y_0_0_yyy, g_x_0_y_0_0_yyz, g_x_0_y_0_0_yz, g_x_0_y_0_0_yzz, g_x_0_y_0_0_zz, g_x_0_y_0_y_xx, g_x_0_y_0_y_xy, g_x_0_y_0_y_xz, g_x_0_y_0_y_yy, g_x_0_y_0_y_yz, g_x_0_y_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_y_xx[k] = -g_x_0_y_0_0_xx[k] * ab_y + g_x_0_y_0_0_xxy[k];

                g_x_0_y_0_y_xy[k] = -g_x_0_y_0_0_xy[k] * ab_y + g_x_0_y_0_0_xyy[k];

                g_x_0_y_0_y_xz[k] = -g_x_0_y_0_0_xz[k] * ab_y + g_x_0_y_0_0_xyz[k];

                g_x_0_y_0_y_yy[k] = -g_x_0_y_0_0_yy[k] * ab_y + g_x_0_y_0_0_yyy[k];

                g_x_0_y_0_y_yz[k] = -g_x_0_y_0_0_yz[k] * ab_y + g_x_0_y_0_0_yyz[k];

                g_x_0_y_0_y_zz[k] = -g_x_0_y_0_0_zz[k] * ab_y + g_x_0_y_0_0_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_z_xx = cbuffer.data(pd_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_y_0_z_xy = cbuffer.data(pd_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_y_0_z_xz = cbuffer.data(pd_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_y_0_z_yy = cbuffer.data(pd_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_y_0_z_yz = cbuffer.data(pd_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_y_0_z_zz = cbuffer.data(pd_geom_1010_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_0_xx, g_x_0_y_0_0_xxz, g_x_0_y_0_0_xy, g_x_0_y_0_0_xyz, g_x_0_y_0_0_xz, g_x_0_y_0_0_xzz, g_x_0_y_0_0_yy, g_x_0_y_0_0_yyz, g_x_0_y_0_0_yz, g_x_0_y_0_0_yzz, g_x_0_y_0_0_zz, g_x_0_y_0_0_zzz, g_x_0_y_0_z_xx, g_x_0_y_0_z_xy, g_x_0_y_0_z_xz, g_x_0_y_0_z_yy, g_x_0_y_0_z_yz, g_x_0_y_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_z_xx[k] = -g_x_0_y_0_0_xx[k] * ab_z + g_x_0_y_0_0_xxz[k];

                g_x_0_y_0_z_xy[k] = -g_x_0_y_0_0_xy[k] * ab_z + g_x_0_y_0_0_xyz[k];

                g_x_0_y_0_z_xz[k] = -g_x_0_y_0_0_xz[k] * ab_z + g_x_0_y_0_0_xzz[k];

                g_x_0_y_0_z_yy[k] = -g_x_0_y_0_0_yy[k] * ab_z + g_x_0_y_0_0_yyz[k];

                g_x_0_y_0_z_yz[k] = -g_x_0_y_0_0_yz[k] * ab_z + g_x_0_y_0_0_yzz[k];

                g_x_0_y_0_z_zz[k] = -g_x_0_y_0_0_zz[k] * ab_z + g_x_0_y_0_0_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_x_xx = cbuffer.data(pd_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_z_0_x_xy = cbuffer.data(pd_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_z_0_x_xz = cbuffer.data(pd_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_z_0_x_yy = cbuffer.data(pd_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_z_0_x_yz = cbuffer.data(pd_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_z_0_x_zz = cbuffer.data(pd_geom_1010_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_xx, g_0_0_z_0_0_xy, g_0_0_z_0_0_xz, g_0_0_z_0_0_yy, g_0_0_z_0_0_yz, g_0_0_z_0_0_zz, g_x_0_z_0_0_xx, g_x_0_z_0_0_xxx, g_x_0_z_0_0_xxy, g_x_0_z_0_0_xxz, g_x_0_z_0_0_xy, g_x_0_z_0_0_xyy, g_x_0_z_0_0_xyz, g_x_0_z_0_0_xz, g_x_0_z_0_0_xzz, g_x_0_z_0_0_yy, g_x_0_z_0_0_yz, g_x_0_z_0_0_zz, g_x_0_z_0_x_xx, g_x_0_z_0_x_xy, g_x_0_z_0_x_xz, g_x_0_z_0_x_yy, g_x_0_z_0_x_yz, g_x_0_z_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_x_xx[k] = -g_0_0_z_0_0_xx[k] - g_x_0_z_0_0_xx[k] * ab_x + g_x_0_z_0_0_xxx[k];

                g_x_0_z_0_x_xy[k] = -g_0_0_z_0_0_xy[k] - g_x_0_z_0_0_xy[k] * ab_x + g_x_0_z_0_0_xxy[k];

                g_x_0_z_0_x_xz[k] = -g_0_0_z_0_0_xz[k] - g_x_0_z_0_0_xz[k] * ab_x + g_x_0_z_0_0_xxz[k];

                g_x_0_z_0_x_yy[k] = -g_0_0_z_0_0_yy[k] - g_x_0_z_0_0_yy[k] * ab_x + g_x_0_z_0_0_xyy[k];

                g_x_0_z_0_x_yz[k] = -g_0_0_z_0_0_yz[k] - g_x_0_z_0_0_yz[k] * ab_x + g_x_0_z_0_0_xyz[k];

                g_x_0_z_0_x_zz[k] = -g_0_0_z_0_0_zz[k] - g_x_0_z_0_0_zz[k] * ab_x + g_x_0_z_0_0_xzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_y_xx = cbuffer.data(pd_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_z_0_y_xy = cbuffer.data(pd_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_z_0_y_xz = cbuffer.data(pd_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_z_0_y_yy = cbuffer.data(pd_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_z_0_y_yz = cbuffer.data(pd_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_z_0_y_zz = cbuffer.data(pd_geom_1010_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_0_xx, g_x_0_z_0_0_xxy, g_x_0_z_0_0_xy, g_x_0_z_0_0_xyy, g_x_0_z_0_0_xyz, g_x_0_z_0_0_xz, g_x_0_z_0_0_yy, g_x_0_z_0_0_yyy, g_x_0_z_0_0_yyz, g_x_0_z_0_0_yz, g_x_0_z_0_0_yzz, g_x_0_z_0_0_zz, g_x_0_z_0_y_xx, g_x_0_z_0_y_xy, g_x_0_z_0_y_xz, g_x_0_z_0_y_yy, g_x_0_z_0_y_yz, g_x_0_z_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_y_xx[k] = -g_x_0_z_0_0_xx[k] * ab_y + g_x_0_z_0_0_xxy[k];

                g_x_0_z_0_y_xy[k] = -g_x_0_z_0_0_xy[k] * ab_y + g_x_0_z_0_0_xyy[k];

                g_x_0_z_0_y_xz[k] = -g_x_0_z_0_0_xz[k] * ab_y + g_x_0_z_0_0_xyz[k];

                g_x_0_z_0_y_yy[k] = -g_x_0_z_0_0_yy[k] * ab_y + g_x_0_z_0_0_yyy[k];

                g_x_0_z_0_y_yz[k] = -g_x_0_z_0_0_yz[k] * ab_y + g_x_0_z_0_0_yyz[k];

                g_x_0_z_0_y_zz[k] = -g_x_0_z_0_0_zz[k] * ab_y + g_x_0_z_0_0_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_z_xx = cbuffer.data(pd_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_z_0_z_xy = cbuffer.data(pd_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_z_0_z_xz = cbuffer.data(pd_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_z_0_z_yy = cbuffer.data(pd_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_z_0_z_yz = cbuffer.data(pd_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_z_0_z_zz = cbuffer.data(pd_geom_1010_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_0_xx, g_x_0_z_0_0_xxz, g_x_0_z_0_0_xy, g_x_0_z_0_0_xyz, g_x_0_z_0_0_xz, g_x_0_z_0_0_xzz, g_x_0_z_0_0_yy, g_x_0_z_0_0_yyz, g_x_0_z_0_0_yz, g_x_0_z_0_0_yzz, g_x_0_z_0_0_zz, g_x_0_z_0_0_zzz, g_x_0_z_0_z_xx, g_x_0_z_0_z_xy, g_x_0_z_0_z_xz, g_x_0_z_0_z_yy, g_x_0_z_0_z_yz, g_x_0_z_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_z_xx[k] = -g_x_0_z_0_0_xx[k] * ab_z + g_x_0_z_0_0_xxz[k];

                g_x_0_z_0_z_xy[k] = -g_x_0_z_0_0_xy[k] * ab_z + g_x_0_z_0_0_xyz[k];

                g_x_0_z_0_z_xz[k] = -g_x_0_z_0_0_xz[k] * ab_z + g_x_0_z_0_0_xzz[k];

                g_x_0_z_0_z_yy[k] = -g_x_0_z_0_0_yy[k] * ab_z + g_x_0_z_0_0_yyz[k];

                g_x_0_z_0_z_yz[k] = -g_x_0_z_0_0_yz[k] * ab_z + g_x_0_z_0_0_yzz[k];

                g_x_0_z_0_z_zz[k] = -g_x_0_z_0_0_zz[k] * ab_z + g_x_0_z_0_0_zzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_x_xx = cbuffer.data(pd_geom_1010_off + 54 * ccomps * dcomps);

            auto g_y_0_x_0_x_xy = cbuffer.data(pd_geom_1010_off + 55 * ccomps * dcomps);

            auto g_y_0_x_0_x_xz = cbuffer.data(pd_geom_1010_off + 56 * ccomps * dcomps);

            auto g_y_0_x_0_x_yy = cbuffer.data(pd_geom_1010_off + 57 * ccomps * dcomps);

            auto g_y_0_x_0_x_yz = cbuffer.data(pd_geom_1010_off + 58 * ccomps * dcomps);

            auto g_y_0_x_0_x_zz = cbuffer.data(pd_geom_1010_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_0_xx, g_y_0_x_0_0_xxx, g_y_0_x_0_0_xxy, g_y_0_x_0_0_xxz, g_y_0_x_0_0_xy, g_y_0_x_0_0_xyy, g_y_0_x_0_0_xyz, g_y_0_x_0_0_xz, g_y_0_x_0_0_xzz, g_y_0_x_0_0_yy, g_y_0_x_0_0_yz, g_y_0_x_0_0_zz, g_y_0_x_0_x_xx, g_y_0_x_0_x_xy, g_y_0_x_0_x_xz, g_y_0_x_0_x_yy, g_y_0_x_0_x_yz, g_y_0_x_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_x_xx[k] = -g_y_0_x_0_0_xx[k] * ab_x + g_y_0_x_0_0_xxx[k];

                g_y_0_x_0_x_xy[k] = -g_y_0_x_0_0_xy[k] * ab_x + g_y_0_x_0_0_xxy[k];

                g_y_0_x_0_x_xz[k] = -g_y_0_x_0_0_xz[k] * ab_x + g_y_0_x_0_0_xxz[k];

                g_y_0_x_0_x_yy[k] = -g_y_0_x_0_0_yy[k] * ab_x + g_y_0_x_0_0_xyy[k];

                g_y_0_x_0_x_yz[k] = -g_y_0_x_0_0_yz[k] * ab_x + g_y_0_x_0_0_xyz[k];

                g_y_0_x_0_x_zz[k] = -g_y_0_x_0_0_zz[k] * ab_x + g_y_0_x_0_0_xzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_y_xx = cbuffer.data(pd_geom_1010_off + 60 * ccomps * dcomps);

            auto g_y_0_x_0_y_xy = cbuffer.data(pd_geom_1010_off + 61 * ccomps * dcomps);

            auto g_y_0_x_0_y_xz = cbuffer.data(pd_geom_1010_off + 62 * ccomps * dcomps);

            auto g_y_0_x_0_y_yy = cbuffer.data(pd_geom_1010_off + 63 * ccomps * dcomps);

            auto g_y_0_x_0_y_yz = cbuffer.data(pd_geom_1010_off + 64 * ccomps * dcomps);

            auto g_y_0_x_0_y_zz = cbuffer.data(pd_geom_1010_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_xx, g_0_0_x_0_0_xy, g_0_0_x_0_0_xz, g_0_0_x_0_0_yy, g_0_0_x_0_0_yz, g_0_0_x_0_0_zz, g_y_0_x_0_0_xx, g_y_0_x_0_0_xxy, g_y_0_x_0_0_xy, g_y_0_x_0_0_xyy, g_y_0_x_0_0_xyz, g_y_0_x_0_0_xz, g_y_0_x_0_0_yy, g_y_0_x_0_0_yyy, g_y_0_x_0_0_yyz, g_y_0_x_0_0_yz, g_y_0_x_0_0_yzz, g_y_0_x_0_0_zz, g_y_0_x_0_y_xx, g_y_0_x_0_y_xy, g_y_0_x_0_y_xz, g_y_0_x_0_y_yy, g_y_0_x_0_y_yz, g_y_0_x_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_y_xx[k] = -g_0_0_x_0_0_xx[k] - g_y_0_x_0_0_xx[k] * ab_y + g_y_0_x_0_0_xxy[k];

                g_y_0_x_0_y_xy[k] = -g_0_0_x_0_0_xy[k] - g_y_0_x_0_0_xy[k] * ab_y + g_y_0_x_0_0_xyy[k];

                g_y_0_x_0_y_xz[k] = -g_0_0_x_0_0_xz[k] - g_y_0_x_0_0_xz[k] * ab_y + g_y_0_x_0_0_xyz[k];

                g_y_0_x_0_y_yy[k] = -g_0_0_x_0_0_yy[k] - g_y_0_x_0_0_yy[k] * ab_y + g_y_0_x_0_0_yyy[k];

                g_y_0_x_0_y_yz[k] = -g_0_0_x_0_0_yz[k] - g_y_0_x_0_0_yz[k] * ab_y + g_y_0_x_0_0_yyz[k];

                g_y_0_x_0_y_zz[k] = -g_0_0_x_0_0_zz[k] - g_y_0_x_0_0_zz[k] * ab_y + g_y_0_x_0_0_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_z_xx = cbuffer.data(pd_geom_1010_off + 66 * ccomps * dcomps);

            auto g_y_0_x_0_z_xy = cbuffer.data(pd_geom_1010_off + 67 * ccomps * dcomps);

            auto g_y_0_x_0_z_xz = cbuffer.data(pd_geom_1010_off + 68 * ccomps * dcomps);

            auto g_y_0_x_0_z_yy = cbuffer.data(pd_geom_1010_off + 69 * ccomps * dcomps);

            auto g_y_0_x_0_z_yz = cbuffer.data(pd_geom_1010_off + 70 * ccomps * dcomps);

            auto g_y_0_x_0_z_zz = cbuffer.data(pd_geom_1010_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_0_xx, g_y_0_x_0_0_xxz, g_y_0_x_0_0_xy, g_y_0_x_0_0_xyz, g_y_0_x_0_0_xz, g_y_0_x_0_0_xzz, g_y_0_x_0_0_yy, g_y_0_x_0_0_yyz, g_y_0_x_0_0_yz, g_y_0_x_0_0_yzz, g_y_0_x_0_0_zz, g_y_0_x_0_0_zzz, g_y_0_x_0_z_xx, g_y_0_x_0_z_xy, g_y_0_x_0_z_xz, g_y_0_x_0_z_yy, g_y_0_x_0_z_yz, g_y_0_x_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_z_xx[k] = -g_y_0_x_0_0_xx[k] * ab_z + g_y_0_x_0_0_xxz[k];

                g_y_0_x_0_z_xy[k] = -g_y_0_x_0_0_xy[k] * ab_z + g_y_0_x_0_0_xyz[k];

                g_y_0_x_0_z_xz[k] = -g_y_0_x_0_0_xz[k] * ab_z + g_y_0_x_0_0_xzz[k];

                g_y_0_x_0_z_yy[k] = -g_y_0_x_0_0_yy[k] * ab_z + g_y_0_x_0_0_yyz[k];

                g_y_0_x_0_z_yz[k] = -g_y_0_x_0_0_yz[k] * ab_z + g_y_0_x_0_0_yzz[k];

                g_y_0_x_0_z_zz[k] = -g_y_0_x_0_0_zz[k] * ab_z + g_y_0_x_0_0_zzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_x_xx = cbuffer.data(pd_geom_1010_off + 72 * ccomps * dcomps);

            auto g_y_0_y_0_x_xy = cbuffer.data(pd_geom_1010_off + 73 * ccomps * dcomps);

            auto g_y_0_y_0_x_xz = cbuffer.data(pd_geom_1010_off + 74 * ccomps * dcomps);

            auto g_y_0_y_0_x_yy = cbuffer.data(pd_geom_1010_off + 75 * ccomps * dcomps);

            auto g_y_0_y_0_x_yz = cbuffer.data(pd_geom_1010_off + 76 * ccomps * dcomps);

            auto g_y_0_y_0_x_zz = cbuffer.data(pd_geom_1010_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_0_xx, g_y_0_y_0_0_xxx, g_y_0_y_0_0_xxy, g_y_0_y_0_0_xxz, g_y_0_y_0_0_xy, g_y_0_y_0_0_xyy, g_y_0_y_0_0_xyz, g_y_0_y_0_0_xz, g_y_0_y_0_0_xzz, g_y_0_y_0_0_yy, g_y_0_y_0_0_yz, g_y_0_y_0_0_zz, g_y_0_y_0_x_xx, g_y_0_y_0_x_xy, g_y_0_y_0_x_xz, g_y_0_y_0_x_yy, g_y_0_y_0_x_yz, g_y_0_y_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_x_xx[k] = -g_y_0_y_0_0_xx[k] * ab_x + g_y_0_y_0_0_xxx[k];

                g_y_0_y_0_x_xy[k] = -g_y_0_y_0_0_xy[k] * ab_x + g_y_0_y_0_0_xxy[k];

                g_y_0_y_0_x_xz[k] = -g_y_0_y_0_0_xz[k] * ab_x + g_y_0_y_0_0_xxz[k];

                g_y_0_y_0_x_yy[k] = -g_y_0_y_0_0_yy[k] * ab_x + g_y_0_y_0_0_xyy[k];

                g_y_0_y_0_x_yz[k] = -g_y_0_y_0_0_yz[k] * ab_x + g_y_0_y_0_0_xyz[k];

                g_y_0_y_0_x_zz[k] = -g_y_0_y_0_0_zz[k] * ab_x + g_y_0_y_0_0_xzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_y_xx = cbuffer.data(pd_geom_1010_off + 78 * ccomps * dcomps);

            auto g_y_0_y_0_y_xy = cbuffer.data(pd_geom_1010_off + 79 * ccomps * dcomps);

            auto g_y_0_y_0_y_xz = cbuffer.data(pd_geom_1010_off + 80 * ccomps * dcomps);

            auto g_y_0_y_0_y_yy = cbuffer.data(pd_geom_1010_off + 81 * ccomps * dcomps);

            auto g_y_0_y_0_y_yz = cbuffer.data(pd_geom_1010_off + 82 * ccomps * dcomps);

            auto g_y_0_y_0_y_zz = cbuffer.data(pd_geom_1010_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_xx, g_0_0_y_0_0_xy, g_0_0_y_0_0_xz, g_0_0_y_0_0_yy, g_0_0_y_0_0_yz, g_0_0_y_0_0_zz, g_y_0_y_0_0_xx, g_y_0_y_0_0_xxy, g_y_0_y_0_0_xy, g_y_0_y_0_0_xyy, g_y_0_y_0_0_xyz, g_y_0_y_0_0_xz, g_y_0_y_0_0_yy, g_y_0_y_0_0_yyy, g_y_0_y_0_0_yyz, g_y_0_y_0_0_yz, g_y_0_y_0_0_yzz, g_y_0_y_0_0_zz, g_y_0_y_0_y_xx, g_y_0_y_0_y_xy, g_y_0_y_0_y_xz, g_y_0_y_0_y_yy, g_y_0_y_0_y_yz, g_y_0_y_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_y_xx[k] = -g_0_0_y_0_0_xx[k] - g_y_0_y_0_0_xx[k] * ab_y + g_y_0_y_0_0_xxy[k];

                g_y_0_y_0_y_xy[k] = -g_0_0_y_0_0_xy[k] - g_y_0_y_0_0_xy[k] * ab_y + g_y_0_y_0_0_xyy[k];

                g_y_0_y_0_y_xz[k] = -g_0_0_y_0_0_xz[k] - g_y_0_y_0_0_xz[k] * ab_y + g_y_0_y_0_0_xyz[k];

                g_y_0_y_0_y_yy[k] = -g_0_0_y_0_0_yy[k] - g_y_0_y_0_0_yy[k] * ab_y + g_y_0_y_0_0_yyy[k];

                g_y_0_y_0_y_yz[k] = -g_0_0_y_0_0_yz[k] - g_y_0_y_0_0_yz[k] * ab_y + g_y_0_y_0_0_yyz[k];

                g_y_0_y_0_y_zz[k] = -g_0_0_y_0_0_zz[k] - g_y_0_y_0_0_zz[k] * ab_y + g_y_0_y_0_0_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_z_xx = cbuffer.data(pd_geom_1010_off + 84 * ccomps * dcomps);

            auto g_y_0_y_0_z_xy = cbuffer.data(pd_geom_1010_off + 85 * ccomps * dcomps);

            auto g_y_0_y_0_z_xz = cbuffer.data(pd_geom_1010_off + 86 * ccomps * dcomps);

            auto g_y_0_y_0_z_yy = cbuffer.data(pd_geom_1010_off + 87 * ccomps * dcomps);

            auto g_y_0_y_0_z_yz = cbuffer.data(pd_geom_1010_off + 88 * ccomps * dcomps);

            auto g_y_0_y_0_z_zz = cbuffer.data(pd_geom_1010_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_0_xx, g_y_0_y_0_0_xxz, g_y_0_y_0_0_xy, g_y_0_y_0_0_xyz, g_y_0_y_0_0_xz, g_y_0_y_0_0_xzz, g_y_0_y_0_0_yy, g_y_0_y_0_0_yyz, g_y_0_y_0_0_yz, g_y_0_y_0_0_yzz, g_y_0_y_0_0_zz, g_y_0_y_0_0_zzz, g_y_0_y_0_z_xx, g_y_0_y_0_z_xy, g_y_0_y_0_z_xz, g_y_0_y_0_z_yy, g_y_0_y_0_z_yz, g_y_0_y_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_z_xx[k] = -g_y_0_y_0_0_xx[k] * ab_z + g_y_0_y_0_0_xxz[k];

                g_y_0_y_0_z_xy[k] = -g_y_0_y_0_0_xy[k] * ab_z + g_y_0_y_0_0_xyz[k];

                g_y_0_y_0_z_xz[k] = -g_y_0_y_0_0_xz[k] * ab_z + g_y_0_y_0_0_xzz[k];

                g_y_0_y_0_z_yy[k] = -g_y_0_y_0_0_yy[k] * ab_z + g_y_0_y_0_0_yyz[k];

                g_y_0_y_0_z_yz[k] = -g_y_0_y_0_0_yz[k] * ab_z + g_y_0_y_0_0_yzz[k];

                g_y_0_y_0_z_zz[k] = -g_y_0_y_0_0_zz[k] * ab_z + g_y_0_y_0_0_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_x_xx = cbuffer.data(pd_geom_1010_off + 90 * ccomps * dcomps);

            auto g_y_0_z_0_x_xy = cbuffer.data(pd_geom_1010_off + 91 * ccomps * dcomps);

            auto g_y_0_z_0_x_xz = cbuffer.data(pd_geom_1010_off + 92 * ccomps * dcomps);

            auto g_y_0_z_0_x_yy = cbuffer.data(pd_geom_1010_off + 93 * ccomps * dcomps);

            auto g_y_0_z_0_x_yz = cbuffer.data(pd_geom_1010_off + 94 * ccomps * dcomps);

            auto g_y_0_z_0_x_zz = cbuffer.data(pd_geom_1010_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_0_xx, g_y_0_z_0_0_xxx, g_y_0_z_0_0_xxy, g_y_0_z_0_0_xxz, g_y_0_z_0_0_xy, g_y_0_z_0_0_xyy, g_y_0_z_0_0_xyz, g_y_0_z_0_0_xz, g_y_0_z_0_0_xzz, g_y_0_z_0_0_yy, g_y_0_z_0_0_yz, g_y_0_z_0_0_zz, g_y_0_z_0_x_xx, g_y_0_z_0_x_xy, g_y_0_z_0_x_xz, g_y_0_z_0_x_yy, g_y_0_z_0_x_yz, g_y_0_z_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_x_xx[k] = -g_y_0_z_0_0_xx[k] * ab_x + g_y_0_z_0_0_xxx[k];

                g_y_0_z_0_x_xy[k] = -g_y_0_z_0_0_xy[k] * ab_x + g_y_0_z_0_0_xxy[k];

                g_y_0_z_0_x_xz[k] = -g_y_0_z_0_0_xz[k] * ab_x + g_y_0_z_0_0_xxz[k];

                g_y_0_z_0_x_yy[k] = -g_y_0_z_0_0_yy[k] * ab_x + g_y_0_z_0_0_xyy[k];

                g_y_0_z_0_x_yz[k] = -g_y_0_z_0_0_yz[k] * ab_x + g_y_0_z_0_0_xyz[k];

                g_y_0_z_0_x_zz[k] = -g_y_0_z_0_0_zz[k] * ab_x + g_y_0_z_0_0_xzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_y_xx = cbuffer.data(pd_geom_1010_off + 96 * ccomps * dcomps);

            auto g_y_0_z_0_y_xy = cbuffer.data(pd_geom_1010_off + 97 * ccomps * dcomps);

            auto g_y_0_z_0_y_xz = cbuffer.data(pd_geom_1010_off + 98 * ccomps * dcomps);

            auto g_y_0_z_0_y_yy = cbuffer.data(pd_geom_1010_off + 99 * ccomps * dcomps);

            auto g_y_0_z_0_y_yz = cbuffer.data(pd_geom_1010_off + 100 * ccomps * dcomps);

            auto g_y_0_z_0_y_zz = cbuffer.data(pd_geom_1010_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_xx, g_0_0_z_0_0_xy, g_0_0_z_0_0_xz, g_0_0_z_0_0_yy, g_0_0_z_0_0_yz, g_0_0_z_0_0_zz, g_y_0_z_0_0_xx, g_y_0_z_0_0_xxy, g_y_0_z_0_0_xy, g_y_0_z_0_0_xyy, g_y_0_z_0_0_xyz, g_y_0_z_0_0_xz, g_y_0_z_0_0_yy, g_y_0_z_0_0_yyy, g_y_0_z_0_0_yyz, g_y_0_z_0_0_yz, g_y_0_z_0_0_yzz, g_y_0_z_0_0_zz, g_y_0_z_0_y_xx, g_y_0_z_0_y_xy, g_y_0_z_0_y_xz, g_y_0_z_0_y_yy, g_y_0_z_0_y_yz, g_y_0_z_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_y_xx[k] = -g_0_0_z_0_0_xx[k] - g_y_0_z_0_0_xx[k] * ab_y + g_y_0_z_0_0_xxy[k];

                g_y_0_z_0_y_xy[k] = -g_0_0_z_0_0_xy[k] - g_y_0_z_0_0_xy[k] * ab_y + g_y_0_z_0_0_xyy[k];

                g_y_0_z_0_y_xz[k] = -g_0_0_z_0_0_xz[k] - g_y_0_z_0_0_xz[k] * ab_y + g_y_0_z_0_0_xyz[k];

                g_y_0_z_0_y_yy[k] = -g_0_0_z_0_0_yy[k] - g_y_0_z_0_0_yy[k] * ab_y + g_y_0_z_0_0_yyy[k];

                g_y_0_z_0_y_yz[k] = -g_0_0_z_0_0_yz[k] - g_y_0_z_0_0_yz[k] * ab_y + g_y_0_z_0_0_yyz[k];

                g_y_0_z_0_y_zz[k] = -g_0_0_z_0_0_zz[k] - g_y_0_z_0_0_zz[k] * ab_y + g_y_0_z_0_0_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_z_xx = cbuffer.data(pd_geom_1010_off + 102 * ccomps * dcomps);

            auto g_y_0_z_0_z_xy = cbuffer.data(pd_geom_1010_off + 103 * ccomps * dcomps);

            auto g_y_0_z_0_z_xz = cbuffer.data(pd_geom_1010_off + 104 * ccomps * dcomps);

            auto g_y_0_z_0_z_yy = cbuffer.data(pd_geom_1010_off + 105 * ccomps * dcomps);

            auto g_y_0_z_0_z_yz = cbuffer.data(pd_geom_1010_off + 106 * ccomps * dcomps);

            auto g_y_0_z_0_z_zz = cbuffer.data(pd_geom_1010_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_0_xx, g_y_0_z_0_0_xxz, g_y_0_z_0_0_xy, g_y_0_z_0_0_xyz, g_y_0_z_0_0_xz, g_y_0_z_0_0_xzz, g_y_0_z_0_0_yy, g_y_0_z_0_0_yyz, g_y_0_z_0_0_yz, g_y_0_z_0_0_yzz, g_y_0_z_0_0_zz, g_y_0_z_0_0_zzz, g_y_0_z_0_z_xx, g_y_0_z_0_z_xy, g_y_0_z_0_z_xz, g_y_0_z_0_z_yy, g_y_0_z_0_z_yz, g_y_0_z_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_z_xx[k] = -g_y_0_z_0_0_xx[k] * ab_z + g_y_0_z_0_0_xxz[k];

                g_y_0_z_0_z_xy[k] = -g_y_0_z_0_0_xy[k] * ab_z + g_y_0_z_0_0_xyz[k];

                g_y_0_z_0_z_xz[k] = -g_y_0_z_0_0_xz[k] * ab_z + g_y_0_z_0_0_xzz[k];

                g_y_0_z_0_z_yy[k] = -g_y_0_z_0_0_yy[k] * ab_z + g_y_0_z_0_0_yyz[k];

                g_y_0_z_0_z_yz[k] = -g_y_0_z_0_0_yz[k] * ab_z + g_y_0_z_0_0_yzz[k];

                g_y_0_z_0_z_zz[k] = -g_y_0_z_0_0_zz[k] * ab_z + g_y_0_z_0_0_zzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_x_xx = cbuffer.data(pd_geom_1010_off + 108 * ccomps * dcomps);

            auto g_z_0_x_0_x_xy = cbuffer.data(pd_geom_1010_off + 109 * ccomps * dcomps);

            auto g_z_0_x_0_x_xz = cbuffer.data(pd_geom_1010_off + 110 * ccomps * dcomps);

            auto g_z_0_x_0_x_yy = cbuffer.data(pd_geom_1010_off + 111 * ccomps * dcomps);

            auto g_z_0_x_0_x_yz = cbuffer.data(pd_geom_1010_off + 112 * ccomps * dcomps);

            auto g_z_0_x_0_x_zz = cbuffer.data(pd_geom_1010_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_0_xx, g_z_0_x_0_0_xxx, g_z_0_x_0_0_xxy, g_z_0_x_0_0_xxz, g_z_0_x_0_0_xy, g_z_0_x_0_0_xyy, g_z_0_x_0_0_xyz, g_z_0_x_0_0_xz, g_z_0_x_0_0_xzz, g_z_0_x_0_0_yy, g_z_0_x_0_0_yz, g_z_0_x_0_0_zz, g_z_0_x_0_x_xx, g_z_0_x_0_x_xy, g_z_0_x_0_x_xz, g_z_0_x_0_x_yy, g_z_0_x_0_x_yz, g_z_0_x_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_x_xx[k] = -g_z_0_x_0_0_xx[k] * ab_x + g_z_0_x_0_0_xxx[k];

                g_z_0_x_0_x_xy[k] = -g_z_0_x_0_0_xy[k] * ab_x + g_z_0_x_0_0_xxy[k];

                g_z_0_x_0_x_xz[k] = -g_z_0_x_0_0_xz[k] * ab_x + g_z_0_x_0_0_xxz[k];

                g_z_0_x_0_x_yy[k] = -g_z_0_x_0_0_yy[k] * ab_x + g_z_0_x_0_0_xyy[k];

                g_z_0_x_0_x_yz[k] = -g_z_0_x_0_0_yz[k] * ab_x + g_z_0_x_0_0_xyz[k];

                g_z_0_x_0_x_zz[k] = -g_z_0_x_0_0_zz[k] * ab_x + g_z_0_x_0_0_xzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_y_xx = cbuffer.data(pd_geom_1010_off + 114 * ccomps * dcomps);

            auto g_z_0_x_0_y_xy = cbuffer.data(pd_geom_1010_off + 115 * ccomps * dcomps);

            auto g_z_0_x_0_y_xz = cbuffer.data(pd_geom_1010_off + 116 * ccomps * dcomps);

            auto g_z_0_x_0_y_yy = cbuffer.data(pd_geom_1010_off + 117 * ccomps * dcomps);

            auto g_z_0_x_0_y_yz = cbuffer.data(pd_geom_1010_off + 118 * ccomps * dcomps);

            auto g_z_0_x_0_y_zz = cbuffer.data(pd_geom_1010_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_0_xx, g_z_0_x_0_0_xxy, g_z_0_x_0_0_xy, g_z_0_x_0_0_xyy, g_z_0_x_0_0_xyz, g_z_0_x_0_0_xz, g_z_0_x_0_0_yy, g_z_0_x_0_0_yyy, g_z_0_x_0_0_yyz, g_z_0_x_0_0_yz, g_z_0_x_0_0_yzz, g_z_0_x_0_0_zz, g_z_0_x_0_y_xx, g_z_0_x_0_y_xy, g_z_0_x_0_y_xz, g_z_0_x_0_y_yy, g_z_0_x_0_y_yz, g_z_0_x_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_y_xx[k] = -g_z_0_x_0_0_xx[k] * ab_y + g_z_0_x_0_0_xxy[k];

                g_z_0_x_0_y_xy[k] = -g_z_0_x_0_0_xy[k] * ab_y + g_z_0_x_0_0_xyy[k];

                g_z_0_x_0_y_xz[k] = -g_z_0_x_0_0_xz[k] * ab_y + g_z_0_x_0_0_xyz[k];

                g_z_0_x_0_y_yy[k] = -g_z_0_x_0_0_yy[k] * ab_y + g_z_0_x_0_0_yyy[k];

                g_z_0_x_0_y_yz[k] = -g_z_0_x_0_0_yz[k] * ab_y + g_z_0_x_0_0_yyz[k];

                g_z_0_x_0_y_zz[k] = -g_z_0_x_0_0_zz[k] * ab_y + g_z_0_x_0_0_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_z_xx = cbuffer.data(pd_geom_1010_off + 120 * ccomps * dcomps);

            auto g_z_0_x_0_z_xy = cbuffer.data(pd_geom_1010_off + 121 * ccomps * dcomps);

            auto g_z_0_x_0_z_xz = cbuffer.data(pd_geom_1010_off + 122 * ccomps * dcomps);

            auto g_z_0_x_0_z_yy = cbuffer.data(pd_geom_1010_off + 123 * ccomps * dcomps);

            auto g_z_0_x_0_z_yz = cbuffer.data(pd_geom_1010_off + 124 * ccomps * dcomps);

            auto g_z_0_x_0_z_zz = cbuffer.data(pd_geom_1010_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_xx, g_0_0_x_0_0_xy, g_0_0_x_0_0_xz, g_0_0_x_0_0_yy, g_0_0_x_0_0_yz, g_0_0_x_0_0_zz, g_z_0_x_0_0_xx, g_z_0_x_0_0_xxz, g_z_0_x_0_0_xy, g_z_0_x_0_0_xyz, g_z_0_x_0_0_xz, g_z_0_x_0_0_xzz, g_z_0_x_0_0_yy, g_z_0_x_0_0_yyz, g_z_0_x_0_0_yz, g_z_0_x_0_0_yzz, g_z_0_x_0_0_zz, g_z_0_x_0_0_zzz, g_z_0_x_0_z_xx, g_z_0_x_0_z_xy, g_z_0_x_0_z_xz, g_z_0_x_0_z_yy, g_z_0_x_0_z_yz, g_z_0_x_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_z_xx[k] = -g_0_0_x_0_0_xx[k] - g_z_0_x_0_0_xx[k] * ab_z + g_z_0_x_0_0_xxz[k];

                g_z_0_x_0_z_xy[k] = -g_0_0_x_0_0_xy[k] - g_z_0_x_0_0_xy[k] * ab_z + g_z_0_x_0_0_xyz[k];

                g_z_0_x_0_z_xz[k] = -g_0_0_x_0_0_xz[k] - g_z_0_x_0_0_xz[k] * ab_z + g_z_0_x_0_0_xzz[k];

                g_z_0_x_0_z_yy[k] = -g_0_0_x_0_0_yy[k] - g_z_0_x_0_0_yy[k] * ab_z + g_z_0_x_0_0_yyz[k];

                g_z_0_x_0_z_yz[k] = -g_0_0_x_0_0_yz[k] - g_z_0_x_0_0_yz[k] * ab_z + g_z_0_x_0_0_yzz[k];

                g_z_0_x_0_z_zz[k] = -g_0_0_x_0_0_zz[k] - g_z_0_x_0_0_zz[k] * ab_z + g_z_0_x_0_0_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_x_xx = cbuffer.data(pd_geom_1010_off + 126 * ccomps * dcomps);

            auto g_z_0_y_0_x_xy = cbuffer.data(pd_geom_1010_off + 127 * ccomps * dcomps);

            auto g_z_0_y_0_x_xz = cbuffer.data(pd_geom_1010_off + 128 * ccomps * dcomps);

            auto g_z_0_y_0_x_yy = cbuffer.data(pd_geom_1010_off + 129 * ccomps * dcomps);

            auto g_z_0_y_0_x_yz = cbuffer.data(pd_geom_1010_off + 130 * ccomps * dcomps);

            auto g_z_0_y_0_x_zz = cbuffer.data(pd_geom_1010_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_0_xx, g_z_0_y_0_0_xxx, g_z_0_y_0_0_xxy, g_z_0_y_0_0_xxz, g_z_0_y_0_0_xy, g_z_0_y_0_0_xyy, g_z_0_y_0_0_xyz, g_z_0_y_0_0_xz, g_z_0_y_0_0_xzz, g_z_0_y_0_0_yy, g_z_0_y_0_0_yz, g_z_0_y_0_0_zz, g_z_0_y_0_x_xx, g_z_0_y_0_x_xy, g_z_0_y_0_x_xz, g_z_0_y_0_x_yy, g_z_0_y_0_x_yz, g_z_0_y_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_x_xx[k] = -g_z_0_y_0_0_xx[k] * ab_x + g_z_0_y_0_0_xxx[k];

                g_z_0_y_0_x_xy[k] = -g_z_0_y_0_0_xy[k] * ab_x + g_z_0_y_0_0_xxy[k];

                g_z_0_y_0_x_xz[k] = -g_z_0_y_0_0_xz[k] * ab_x + g_z_0_y_0_0_xxz[k];

                g_z_0_y_0_x_yy[k] = -g_z_0_y_0_0_yy[k] * ab_x + g_z_0_y_0_0_xyy[k];

                g_z_0_y_0_x_yz[k] = -g_z_0_y_0_0_yz[k] * ab_x + g_z_0_y_0_0_xyz[k];

                g_z_0_y_0_x_zz[k] = -g_z_0_y_0_0_zz[k] * ab_x + g_z_0_y_0_0_xzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_y_xx = cbuffer.data(pd_geom_1010_off + 132 * ccomps * dcomps);

            auto g_z_0_y_0_y_xy = cbuffer.data(pd_geom_1010_off + 133 * ccomps * dcomps);

            auto g_z_0_y_0_y_xz = cbuffer.data(pd_geom_1010_off + 134 * ccomps * dcomps);

            auto g_z_0_y_0_y_yy = cbuffer.data(pd_geom_1010_off + 135 * ccomps * dcomps);

            auto g_z_0_y_0_y_yz = cbuffer.data(pd_geom_1010_off + 136 * ccomps * dcomps);

            auto g_z_0_y_0_y_zz = cbuffer.data(pd_geom_1010_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_0_xx, g_z_0_y_0_0_xxy, g_z_0_y_0_0_xy, g_z_0_y_0_0_xyy, g_z_0_y_0_0_xyz, g_z_0_y_0_0_xz, g_z_0_y_0_0_yy, g_z_0_y_0_0_yyy, g_z_0_y_0_0_yyz, g_z_0_y_0_0_yz, g_z_0_y_0_0_yzz, g_z_0_y_0_0_zz, g_z_0_y_0_y_xx, g_z_0_y_0_y_xy, g_z_0_y_0_y_xz, g_z_0_y_0_y_yy, g_z_0_y_0_y_yz, g_z_0_y_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_y_xx[k] = -g_z_0_y_0_0_xx[k] * ab_y + g_z_0_y_0_0_xxy[k];

                g_z_0_y_0_y_xy[k] = -g_z_0_y_0_0_xy[k] * ab_y + g_z_0_y_0_0_xyy[k];

                g_z_0_y_0_y_xz[k] = -g_z_0_y_0_0_xz[k] * ab_y + g_z_0_y_0_0_xyz[k];

                g_z_0_y_0_y_yy[k] = -g_z_0_y_0_0_yy[k] * ab_y + g_z_0_y_0_0_yyy[k];

                g_z_0_y_0_y_yz[k] = -g_z_0_y_0_0_yz[k] * ab_y + g_z_0_y_0_0_yyz[k];

                g_z_0_y_0_y_zz[k] = -g_z_0_y_0_0_zz[k] * ab_y + g_z_0_y_0_0_yzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_z_xx = cbuffer.data(pd_geom_1010_off + 138 * ccomps * dcomps);

            auto g_z_0_y_0_z_xy = cbuffer.data(pd_geom_1010_off + 139 * ccomps * dcomps);

            auto g_z_0_y_0_z_xz = cbuffer.data(pd_geom_1010_off + 140 * ccomps * dcomps);

            auto g_z_0_y_0_z_yy = cbuffer.data(pd_geom_1010_off + 141 * ccomps * dcomps);

            auto g_z_0_y_0_z_yz = cbuffer.data(pd_geom_1010_off + 142 * ccomps * dcomps);

            auto g_z_0_y_0_z_zz = cbuffer.data(pd_geom_1010_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_xx, g_0_0_y_0_0_xy, g_0_0_y_0_0_xz, g_0_0_y_0_0_yy, g_0_0_y_0_0_yz, g_0_0_y_0_0_zz, g_z_0_y_0_0_xx, g_z_0_y_0_0_xxz, g_z_0_y_0_0_xy, g_z_0_y_0_0_xyz, g_z_0_y_0_0_xz, g_z_0_y_0_0_xzz, g_z_0_y_0_0_yy, g_z_0_y_0_0_yyz, g_z_0_y_0_0_yz, g_z_0_y_0_0_yzz, g_z_0_y_0_0_zz, g_z_0_y_0_0_zzz, g_z_0_y_0_z_xx, g_z_0_y_0_z_xy, g_z_0_y_0_z_xz, g_z_0_y_0_z_yy, g_z_0_y_0_z_yz, g_z_0_y_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_z_xx[k] = -g_0_0_y_0_0_xx[k] - g_z_0_y_0_0_xx[k] * ab_z + g_z_0_y_0_0_xxz[k];

                g_z_0_y_0_z_xy[k] = -g_0_0_y_0_0_xy[k] - g_z_0_y_0_0_xy[k] * ab_z + g_z_0_y_0_0_xyz[k];

                g_z_0_y_0_z_xz[k] = -g_0_0_y_0_0_xz[k] - g_z_0_y_0_0_xz[k] * ab_z + g_z_0_y_0_0_xzz[k];

                g_z_0_y_0_z_yy[k] = -g_0_0_y_0_0_yy[k] - g_z_0_y_0_0_yy[k] * ab_z + g_z_0_y_0_0_yyz[k];

                g_z_0_y_0_z_yz[k] = -g_0_0_y_0_0_yz[k] - g_z_0_y_0_0_yz[k] * ab_z + g_z_0_y_0_0_yzz[k];

                g_z_0_y_0_z_zz[k] = -g_0_0_y_0_0_zz[k] - g_z_0_y_0_0_zz[k] * ab_z + g_z_0_y_0_0_zzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_x_xx = cbuffer.data(pd_geom_1010_off + 144 * ccomps * dcomps);

            auto g_z_0_z_0_x_xy = cbuffer.data(pd_geom_1010_off + 145 * ccomps * dcomps);

            auto g_z_0_z_0_x_xz = cbuffer.data(pd_geom_1010_off + 146 * ccomps * dcomps);

            auto g_z_0_z_0_x_yy = cbuffer.data(pd_geom_1010_off + 147 * ccomps * dcomps);

            auto g_z_0_z_0_x_yz = cbuffer.data(pd_geom_1010_off + 148 * ccomps * dcomps);

            auto g_z_0_z_0_x_zz = cbuffer.data(pd_geom_1010_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_0_xx, g_z_0_z_0_0_xxx, g_z_0_z_0_0_xxy, g_z_0_z_0_0_xxz, g_z_0_z_0_0_xy, g_z_0_z_0_0_xyy, g_z_0_z_0_0_xyz, g_z_0_z_0_0_xz, g_z_0_z_0_0_xzz, g_z_0_z_0_0_yy, g_z_0_z_0_0_yz, g_z_0_z_0_0_zz, g_z_0_z_0_x_xx, g_z_0_z_0_x_xy, g_z_0_z_0_x_xz, g_z_0_z_0_x_yy, g_z_0_z_0_x_yz, g_z_0_z_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_x_xx[k] = -g_z_0_z_0_0_xx[k] * ab_x + g_z_0_z_0_0_xxx[k];

                g_z_0_z_0_x_xy[k] = -g_z_0_z_0_0_xy[k] * ab_x + g_z_0_z_0_0_xxy[k];

                g_z_0_z_0_x_xz[k] = -g_z_0_z_0_0_xz[k] * ab_x + g_z_0_z_0_0_xxz[k];

                g_z_0_z_0_x_yy[k] = -g_z_0_z_0_0_yy[k] * ab_x + g_z_0_z_0_0_xyy[k];

                g_z_0_z_0_x_yz[k] = -g_z_0_z_0_0_yz[k] * ab_x + g_z_0_z_0_0_xyz[k];

                g_z_0_z_0_x_zz[k] = -g_z_0_z_0_0_zz[k] * ab_x + g_z_0_z_0_0_xzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_y_xx = cbuffer.data(pd_geom_1010_off + 150 * ccomps * dcomps);

            auto g_z_0_z_0_y_xy = cbuffer.data(pd_geom_1010_off + 151 * ccomps * dcomps);

            auto g_z_0_z_0_y_xz = cbuffer.data(pd_geom_1010_off + 152 * ccomps * dcomps);

            auto g_z_0_z_0_y_yy = cbuffer.data(pd_geom_1010_off + 153 * ccomps * dcomps);

            auto g_z_0_z_0_y_yz = cbuffer.data(pd_geom_1010_off + 154 * ccomps * dcomps);

            auto g_z_0_z_0_y_zz = cbuffer.data(pd_geom_1010_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_0_xx, g_z_0_z_0_0_xxy, g_z_0_z_0_0_xy, g_z_0_z_0_0_xyy, g_z_0_z_0_0_xyz, g_z_0_z_0_0_xz, g_z_0_z_0_0_yy, g_z_0_z_0_0_yyy, g_z_0_z_0_0_yyz, g_z_0_z_0_0_yz, g_z_0_z_0_0_yzz, g_z_0_z_0_0_zz, g_z_0_z_0_y_xx, g_z_0_z_0_y_xy, g_z_0_z_0_y_xz, g_z_0_z_0_y_yy, g_z_0_z_0_y_yz, g_z_0_z_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_y_xx[k] = -g_z_0_z_0_0_xx[k] * ab_y + g_z_0_z_0_0_xxy[k];

                g_z_0_z_0_y_xy[k] = -g_z_0_z_0_0_xy[k] * ab_y + g_z_0_z_0_0_xyy[k];

                g_z_0_z_0_y_xz[k] = -g_z_0_z_0_0_xz[k] * ab_y + g_z_0_z_0_0_xyz[k];

                g_z_0_z_0_y_yy[k] = -g_z_0_z_0_0_yy[k] * ab_y + g_z_0_z_0_0_yyy[k];

                g_z_0_z_0_y_yz[k] = -g_z_0_z_0_0_yz[k] * ab_y + g_z_0_z_0_0_yyz[k];

                g_z_0_z_0_y_zz[k] = -g_z_0_z_0_0_zz[k] * ab_y + g_z_0_z_0_0_yzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_z_xx = cbuffer.data(pd_geom_1010_off + 156 * ccomps * dcomps);

            auto g_z_0_z_0_z_xy = cbuffer.data(pd_geom_1010_off + 157 * ccomps * dcomps);

            auto g_z_0_z_0_z_xz = cbuffer.data(pd_geom_1010_off + 158 * ccomps * dcomps);

            auto g_z_0_z_0_z_yy = cbuffer.data(pd_geom_1010_off + 159 * ccomps * dcomps);

            auto g_z_0_z_0_z_yz = cbuffer.data(pd_geom_1010_off + 160 * ccomps * dcomps);

            auto g_z_0_z_0_z_zz = cbuffer.data(pd_geom_1010_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_xx, g_0_0_z_0_0_xy, g_0_0_z_0_0_xz, g_0_0_z_0_0_yy, g_0_0_z_0_0_yz, g_0_0_z_0_0_zz, g_z_0_z_0_0_xx, g_z_0_z_0_0_xxz, g_z_0_z_0_0_xy, g_z_0_z_0_0_xyz, g_z_0_z_0_0_xz, g_z_0_z_0_0_xzz, g_z_0_z_0_0_yy, g_z_0_z_0_0_yyz, g_z_0_z_0_0_yz, g_z_0_z_0_0_yzz, g_z_0_z_0_0_zz, g_z_0_z_0_0_zzz, g_z_0_z_0_z_xx, g_z_0_z_0_z_xy, g_z_0_z_0_z_xz, g_z_0_z_0_z_yy, g_z_0_z_0_z_yz, g_z_0_z_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_z_xx[k] = -g_0_0_z_0_0_xx[k] - g_z_0_z_0_0_xx[k] * ab_z + g_z_0_z_0_0_xxz[k];

                g_z_0_z_0_z_xy[k] = -g_0_0_z_0_0_xy[k] - g_z_0_z_0_0_xy[k] * ab_z + g_z_0_z_0_0_xyz[k];

                g_z_0_z_0_z_xz[k] = -g_0_0_z_0_0_xz[k] - g_z_0_z_0_0_xz[k] * ab_z + g_z_0_z_0_0_xzz[k];

                g_z_0_z_0_z_yy[k] = -g_0_0_z_0_0_yy[k] - g_z_0_z_0_0_yy[k] * ab_z + g_z_0_z_0_0_yyz[k];

                g_z_0_z_0_z_yz[k] = -g_0_0_z_0_0_yz[k] - g_z_0_z_0_0_yz[k] * ab_z + g_z_0_z_0_0_yzz[k];

                g_z_0_z_0_z_zz[k] = -g_0_0_z_0_0_zz[k] - g_z_0_z_0_0_zz[k] * ab_z + g_z_0_z_0_0_zzz[k];
            }
        }
    }
}

} // erirec namespace

