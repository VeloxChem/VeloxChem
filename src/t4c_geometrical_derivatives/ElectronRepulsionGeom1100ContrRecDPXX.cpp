#include "ElectronRepulsionGeom1100ContrRecDPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_dpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_dpxx,
                                            const size_t idx_geom_01_ppxx,
                                            const size_t idx_geom_10_ppxx,
                                            const size_t idx_geom_11_ppxx,
                                            const size_t idx_geom_11_pdxx,
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
            /// Set up components of auxilary buffer : PPSS

            const auto pp_geom_01_off = idx_geom_01_ppxx + i * dcomps + j;

            auto g_0_x_x_x = cbuffer.data(pp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_y = cbuffer.data(pp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_z = cbuffer.data(pp_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_y_x = cbuffer.data(pp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_y_y = cbuffer.data(pp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_y_z = cbuffer.data(pp_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_z_x = cbuffer.data(pp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_z_y = cbuffer.data(pp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_z_z = cbuffer.data(pp_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_y_x_x = cbuffer.data(pp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_y_x_y = cbuffer.data(pp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_y_x_z = cbuffer.data(pp_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_y_y_x = cbuffer.data(pp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_y_y_y = cbuffer.data(pp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_y_y_z = cbuffer.data(pp_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_y_z_x = cbuffer.data(pp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_y_z_y = cbuffer.data(pp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_y_z_z = cbuffer.data(pp_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_z_x_x = cbuffer.data(pp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_z_x_y = cbuffer.data(pp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_z_x_z = cbuffer.data(pp_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_z_y_x = cbuffer.data(pp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_z_y_y = cbuffer.data(pp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_z_y_z = cbuffer.data(pp_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_z_z_x = cbuffer.data(pp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_z_z_y = cbuffer.data(pp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_z_z_z = cbuffer.data(pp_geom_01_off + 26 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PPSS

            const auto pp_geom_10_off = idx_geom_10_ppxx + i * dcomps + j;

            auto g_x_0_x_x = cbuffer.data(pp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_y = cbuffer.data(pp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_z = cbuffer.data(pp_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_y_x = cbuffer.data(pp_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_y_y = cbuffer.data(pp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_y_z = cbuffer.data(pp_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_z_x = cbuffer.data(pp_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_z_y = cbuffer.data(pp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_z_z = cbuffer.data(pp_geom_10_off + 8 * ccomps * dcomps);

            auto g_y_0_x_x = cbuffer.data(pp_geom_10_off + 9 * ccomps * dcomps);

            auto g_y_0_x_y = cbuffer.data(pp_geom_10_off + 10 * ccomps * dcomps);

            auto g_y_0_x_z = cbuffer.data(pp_geom_10_off + 11 * ccomps * dcomps);

            auto g_y_0_y_x = cbuffer.data(pp_geom_10_off + 12 * ccomps * dcomps);

            auto g_y_0_y_y = cbuffer.data(pp_geom_10_off + 13 * ccomps * dcomps);

            auto g_y_0_y_z = cbuffer.data(pp_geom_10_off + 14 * ccomps * dcomps);

            auto g_y_0_z_x = cbuffer.data(pp_geom_10_off + 15 * ccomps * dcomps);

            auto g_y_0_z_y = cbuffer.data(pp_geom_10_off + 16 * ccomps * dcomps);

            auto g_y_0_z_z = cbuffer.data(pp_geom_10_off + 17 * ccomps * dcomps);

            auto g_z_0_x_x = cbuffer.data(pp_geom_10_off + 18 * ccomps * dcomps);

            auto g_z_0_x_y = cbuffer.data(pp_geom_10_off + 19 * ccomps * dcomps);

            auto g_z_0_x_z = cbuffer.data(pp_geom_10_off + 20 * ccomps * dcomps);

            auto g_z_0_y_x = cbuffer.data(pp_geom_10_off + 21 * ccomps * dcomps);

            auto g_z_0_y_y = cbuffer.data(pp_geom_10_off + 22 * ccomps * dcomps);

            auto g_z_0_y_z = cbuffer.data(pp_geom_10_off + 23 * ccomps * dcomps);

            auto g_z_0_z_x = cbuffer.data(pp_geom_10_off + 24 * ccomps * dcomps);

            auto g_z_0_z_y = cbuffer.data(pp_geom_10_off + 25 * ccomps * dcomps);

            auto g_z_0_z_z = cbuffer.data(pp_geom_10_off + 26 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PPSS

            const auto pp_geom_11_off = idx_geom_11_ppxx + i * dcomps + j;

            auto g_x_x_x_x = cbuffer.data(pp_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_y = cbuffer.data(pp_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_z = cbuffer.data(pp_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_y_x = cbuffer.data(pp_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_y_y = cbuffer.data(pp_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_y_z = cbuffer.data(pp_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_z_x = cbuffer.data(pp_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_z_y = cbuffer.data(pp_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_z_z = cbuffer.data(pp_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_y_x_x = cbuffer.data(pp_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_y_x_y = cbuffer.data(pp_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_y_x_z = cbuffer.data(pp_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_y_y_x = cbuffer.data(pp_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_y_y_y = cbuffer.data(pp_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_y_y_z = cbuffer.data(pp_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_y_z_x = cbuffer.data(pp_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_y_z_y = cbuffer.data(pp_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_y_z_z = cbuffer.data(pp_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_z_x_x = cbuffer.data(pp_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_z_x_y = cbuffer.data(pp_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_z_x_z = cbuffer.data(pp_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_z_y_x = cbuffer.data(pp_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_z_y_y = cbuffer.data(pp_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_z_y_z = cbuffer.data(pp_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_z_z_x = cbuffer.data(pp_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_z_z_y = cbuffer.data(pp_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_z_z_z = cbuffer.data(pp_geom_11_off + 26 * ccomps * dcomps);

            auto g_y_x_x_x = cbuffer.data(pp_geom_11_off + 27 * ccomps * dcomps);

            auto g_y_x_x_y = cbuffer.data(pp_geom_11_off + 28 * ccomps * dcomps);

            auto g_y_x_x_z = cbuffer.data(pp_geom_11_off + 29 * ccomps * dcomps);

            auto g_y_x_y_x = cbuffer.data(pp_geom_11_off + 30 * ccomps * dcomps);

            auto g_y_x_y_y = cbuffer.data(pp_geom_11_off + 31 * ccomps * dcomps);

            auto g_y_x_y_z = cbuffer.data(pp_geom_11_off + 32 * ccomps * dcomps);

            auto g_y_x_z_x = cbuffer.data(pp_geom_11_off + 33 * ccomps * dcomps);

            auto g_y_x_z_y = cbuffer.data(pp_geom_11_off + 34 * ccomps * dcomps);

            auto g_y_x_z_z = cbuffer.data(pp_geom_11_off + 35 * ccomps * dcomps);

            auto g_y_y_x_x = cbuffer.data(pp_geom_11_off + 36 * ccomps * dcomps);

            auto g_y_y_x_y = cbuffer.data(pp_geom_11_off + 37 * ccomps * dcomps);

            auto g_y_y_x_z = cbuffer.data(pp_geom_11_off + 38 * ccomps * dcomps);

            auto g_y_y_y_x = cbuffer.data(pp_geom_11_off + 39 * ccomps * dcomps);

            auto g_y_y_y_y = cbuffer.data(pp_geom_11_off + 40 * ccomps * dcomps);

            auto g_y_y_y_z = cbuffer.data(pp_geom_11_off + 41 * ccomps * dcomps);

            auto g_y_y_z_x = cbuffer.data(pp_geom_11_off + 42 * ccomps * dcomps);

            auto g_y_y_z_y = cbuffer.data(pp_geom_11_off + 43 * ccomps * dcomps);

            auto g_y_y_z_z = cbuffer.data(pp_geom_11_off + 44 * ccomps * dcomps);

            auto g_y_z_x_x = cbuffer.data(pp_geom_11_off + 45 * ccomps * dcomps);

            auto g_y_z_x_y = cbuffer.data(pp_geom_11_off + 46 * ccomps * dcomps);

            auto g_y_z_x_z = cbuffer.data(pp_geom_11_off + 47 * ccomps * dcomps);

            auto g_y_z_y_x = cbuffer.data(pp_geom_11_off + 48 * ccomps * dcomps);

            auto g_y_z_y_y = cbuffer.data(pp_geom_11_off + 49 * ccomps * dcomps);

            auto g_y_z_y_z = cbuffer.data(pp_geom_11_off + 50 * ccomps * dcomps);

            auto g_y_z_z_x = cbuffer.data(pp_geom_11_off + 51 * ccomps * dcomps);

            auto g_y_z_z_y = cbuffer.data(pp_geom_11_off + 52 * ccomps * dcomps);

            auto g_y_z_z_z = cbuffer.data(pp_geom_11_off + 53 * ccomps * dcomps);

            auto g_z_x_x_x = cbuffer.data(pp_geom_11_off + 54 * ccomps * dcomps);

            auto g_z_x_x_y = cbuffer.data(pp_geom_11_off + 55 * ccomps * dcomps);

            auto g_z_x_x_z = cbuffer.data(pp_geom_11_off + 56 * ccomps * dcomps);

            auto g_z_x_y_x = cbuffer.data(pp_geom_11_off + 57 * ccomps * dcomps);

            auto g_z_x_y_y = cbuffer.data(pp_geom_11_off + 58 * ccomps * dcomps);

            auto g_z_x_y_z = cbuffer.data(pp_geom_11_off + 59 * ccomps * dcomps);

            auto g_z_x_z_x = cbuffer.data(pp_geom_11_off + 60 * ccomps * dcomps);

            auto g_z_x_z_y = cbuffer.data(pp_geom_11_off + 61 * ccomps * dcomps);

            auto g_z_x_z_z = cbuffer.data(pp_geom_11_off + 62 * ccomps * dcomps);

            auto g_z_y_x_x = cbuffer.data(pp_geom_11_off + 63 * ccomps * dcomps);

            auto g_z_y_x_y = cbuffer.data(pp_geom_11_off + 64 * ccomps * dcomps);

            auto g_z_y_x_z = cbuffer.data(pp_geom_11_off + 65 * ccomps * dcomps);

            auto g_z_y_y_x = cbuffer.data(pp_geom_11_off + 66 * ccomps * dcomps);

            auto g_z_y_y_y = cbuffer.data(pp_geom_11_off + 67 * ccomps * dcomps);

            auto g_z_y_y_z = cbuffer.data(pp_geom_11_off + 68 * ccomps * dcomps);

            auto g_z_y_z_x = cbuffer.data(pp_geom_11_off + 69 * ccomps * dcomps);

            auto g_z_y_z_y = cbuffer.data(pp_geom_11_off + 70 * ccomps * dcomps);

            auto g_z_y_z_z = cbuffer.data(pp_geom_11_off + 71 * ccomps * dcomps);

            auto g_z_z_x_x = cbuffer.data(pp_geom_11_off + 72 * ccomps * dcomps);

            auto g_z_z_x_y = cbuffer.data(pp_geom_11_off + 73 * ccomps * dcomps);

            auto g_z_z_x_z = cbuffer.data(pp_geom_11_off + 74 * ccomps * dcomps);

            auto g_z_z_y_x = cbuffer.data(pp_geom_11_off + 75 * ccomps * dcomps);

            auto g_z_z_y_y = cbuffer.data(pp_geom_11_off + 76 * ccomps * dcomps);

            auto g_z_z_y_z = cbuffer.data(pp_geom_11_off + 77 * ccomps * dcomps);

            auto g_z_z_z_x = cbuffer.data(pp_geom_11_off + 78 * ccomps * dcomps);

            auto g_z_z_z_y = cbuffer.data(pp_geom_11_off + 79 * ccomps * dcomps);

            auto g_z_z_z_z = cbuffer.data(pp_geom_11_off + 80 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PDSS

            const auto pd_geom_11_off = idx_geom_11_pdxx + i * dcomps + j;

            auto g_x_x_x_xx = cbuffer.data(pd_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xy = cbuffer.data(pd_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xz = cbuffer.data(pd_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_yy = cbuffer.data(pd_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_yz = cbuffer.data(pd_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_zz = cbuffer.data(pd_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_y_xy = cbuffer.data(pd_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_y_yy = cbuffer.data(pd_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_y_yz = cbuffer.data(pd_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_z_xy = cbuffer.data(pd_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_z_xz = cbuffer.data(pd_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_z_yy = cbuffer.data(pd_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_z_yz = cbuffer.data(pd_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_z_zz = cbuffer.data(pd_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_y_x_xx = cbuffer.data(pd_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_y_x_xy = cbuffer.data(pd_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_y_x_xz = cbuffer.data(pd_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_y_x_yz = cbuffer.data(pd_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_y_x_zz = cbuffer.data(pd_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_y_y_xx = cbuffer.data(pd_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_y_y_xy = cbuffer.data(pd_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_y_y_xz = cbuffer.data(pd_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_y_y_yy = cbuffer.data(pd_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_y_y_yz = cbuffer.data(pd_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_y_y_zz = cbuffer.data(pd_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_y_z_xz = cbuffer.data(pd_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_y_z_yz = cbuffer.data(pd_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_z_zz = cbuffer.data(pd_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_z_x_xx = cbuffer.data(pd_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_z_x_xy = cbuffer.data(pd_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_z_x_xz = cbuffer.data(pd_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_z_x_yy = cbuffer.data(pd_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_z_x_yz = cbuffer.data(pd_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_z_y_xy = cbuffer.data(pd_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_z_y_yy = cbuffer.data(pd_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_z_y_yz = cbuffer.data(pd_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_z_z_xx = cbuffer.data(pd_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_z_z_xy = cbuffer.data(pd_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_z_z_xz = cbuffer.data(pd_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_z_z_yy = cbuffer.data(pd_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_z_z_yz = cbuffer.data(pd_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_z_z_zz = cbuffer.data(pd_geom_11_off + 53 * ccomps * dcomps);

            auto g_y_x_x_xx = cbuffer.data(pd_geom_11_off + 54 * ccomps * dcomps);

            auto g_y_x_x_xy = cbuffer.data(pd_geom_11_off + 55 * ccomps * dcomps);

            auto g_y_x_x_xz = cbuffer.data(pd_geom_11_off + 56 * ccomps * dcomps);

            auto g_y_x_x_yz = cbuffer.data(pd_geom_11_off + 58 * ccomps * dcomps);

            auto g_y_x_x_zz = cbuffer.data(pd_geom_11_off + 59 * ccomps * dcomps);

            auto g_y_x_y_xx = cbuffer.data(pd_geom_11_off + 60 * ccomps * dcomps);

            auto g_y_x_y_xy = cbuffer.data(pd_geom_11_off + 61 * ccomps * dcomps);

            auto g_y_x_y_xz = cbuffer.data(pd_geom_11_off + 62 * ccomps * dcomps);

            auto g_y_x_y_yy = cbuffer.data(pd_geom_11_off + 63 * ccomps * dcomps);

            auto g_y_x_y_yz = cbuffer.data(pd_geom_11_off + 64 * ccomps * dcomps);

            auto g_y_x_y_zz = cbuffer.data(pd_geom_11_off + 65 * ccomps * dcomps);

            auto g_y_x_z_xz = cbuffer.data(pd_geom_11_off + 68 * ccomps * dcomps);

            auto g_y_x_z_yz = cbuffer.data(pd_geom_11_off + 70 * ccomps * dcomps);

            auto g_y_x_z_zz = cbuffer.data(pd_geom_11_off + 71 * ccomps * dcomps);

            auto g_y_y_x_xx = cbuffer.data(pd_geom_11_off + 72 * ccomps * dcomps);

            auto g_y_y_x_xy = cbuffer.data(pd_geom_11_off + 73 * ccomps * dcomps);

            auto g_y_y_x_xz = cbuffer.data(pd_geom_11_off + 74 * ccomps * dcomps);

            auto g_y_y_y_xx = cbuffer.data(pd_geom_11_off + 78 * ccomps * dcomps);

            auto g_y_y_y_xy = cbuffer.data(pd_geom_11_off + 79 * ccomps * dcomps);

            auto g_y_y_y_xz = cbuffer.data(pd_geom_11_off + 80 * ccomps * dcomps);

            auto g_y_y_y_yy = cbuffer.data(pd_geom_11_off + 81 * ccomps * dcomps);

            auto g_y_y_y_yz = cbuffer.data(pd_geom_11_off + 82 * ccomps * dcomps);

            auto g_y_y_y_zz = cbuffer.data(pd_geom_11_off + 83 * ccomps * dcomps);

            auto g_y_y_z_xx = cbuffer.data(pd_geom_11_off + 84 * ccomps * dcomps);

            auto g_y_y_z_xy = cbuffer.data(pd_geom_11_off + 85 * ccomps * dcomps);

            auto g_y_y_z_xz = cbuffer.data(pd_geom_11_off + 86 * ccomps * dcomps);

            auto g_y_y_z_yz = cbuffer.data(pd_geom_11_off + 88 * ccomps * dcomps);

            auto g_y_y_z_zz = cbuffer.data(pd_geom_11_off + 89 * ccomps * dcomps);

            auto g_y_z_x_xx = cbuffer.data(pd_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_z_x_xy = cbuffer.data(pd_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_z_x_xz = cbuffer.data(pd_geom_11_off + 92 * ccomps * dcomps);

            auto g_y_z_y_xx = cbuffer.data(pd_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_z_y_xy = cbuffer.data(pd_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_z_y_xz = cbuffer.data(pd_geom_11_off + 98 * ccomps * dcomps);

            auto g_y_z_y_yy = cbuffer.data(pd_geom_11_off + 99 * ccomps * dcomps);

            auto g_y_z_y_yz = cbuffer.data(pd_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_z_z_xx = cbuffer.data(pd_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_z_z_xy = cbuffer.data(pd_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_z_z_xz = cbuffer.data(pd_geom_11_off + 104 * ccomps * dcomps);

            auto g_y_z_z_yy = cbuffer.data(pd_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_z_z_yz = cbuffer.data(pd_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_z_z_zz = cbuffer.data(pd_geom_11_off + 107 * ccomps * dcomps);

            auto g_z_x_x_xx = cbuffer.data(pd_geom_11_off + 108 * ccomps * dcomps);

            auto g_z_x_x_xy = cbuffer.data(pd_geom_11_off + 109 * ccomps * dcomps);

            auto g_z_x_x_xz = cbuffer.data(pd_geom_11_off + 110 * ccomps * dcomps);

            auto g_z_x_x_yy = cbuffer.data(pd_geom_11_off + 111 * ccomps * dcomps);

            auto g_z_x_x_yz = cbuffer.data(pd_geom_11_off + 112 * ccomps * dcomps);

            auto g_z_x_y_xy = cbuffer.data(pd_geom_11_off + 115 * ccomps * dcomps);

            auto g_z_x_y_yy = cbuffer.data(pd_geom_11_off + 117 * ccomps * dcomps);

            auto g_z_x_y_yz = cbuffer.data(pd_geom_11_off + 118 * ccomps * dcomps);

            auto g_z_x_z_xx = cbuffer.data(pd_geom_11_off + 120 * ccomps * dcomps);

            auto g_z_x_z_xy = cbuffer.data(pd_geom_11_off + 121 * ccomps * dcomps);

            auto g_z_x_z_xz = cbuffer.data(pd_geom_11_off + 122 * ccomps * dcomps);

            auto g_z_x_z_yy = cbuffer.data(pd_geom_11_off + 123 * ccomps * dcomps);

            auto g_z_x_z_yz = cbuffer.data(pd_geom_11_off + 124 * ccomps * dcomps);

            auto g_z_x_z_zz = cbuffer.data(pd_geom_11_off + 125 * ccomps * dcomps);

            auto g_z_y_x_xx = cbuffer.data(pd_geom_11_off + 126 * ccomps * dcomps);

            auto g_z_y_x_xy = cbuffer.data(pd_geom_11_off + 127 * ccomps * dcomps);

            auto g_z_y_x_xz = cbuffer.data(pd_geom_11_off + 128 * ccomps * dcomps);

            auto g_z_y_y_xx = cbuffer.data(pd_geom_11_off + 132 * ccomps * dcomps);

            auto g_z_y_y_xy = cbuffer.data(pd_geom_11_off + 133 * ccomps * dcomps);

            auto g_z_y_y_xz = cbuffer.data(pd_geom_11_off + 134 * ccomps * dcomps);

            auto g_z_y_y_yy = cbuffer.data(pd_geom_11_off + 135 * ccomps * dcomps);

            auto g_z_y_y_yz = cbuffer.data(pd_geom_11_off + 136 * ccomps * dcomps);

            auto g_z_y_z_xx = cbuffer.data(pd_geom_11_off + 138 * ccomps * dcomps);

            auto g_z_y_z_xy = cbuffer.data(pd_geom_11_off + 139 * ccomps * dcomps);

            auto g_z_y_z_xz = cbuffer.data(pd_geom_11_off + 140 * ccomps * dcomps);

            auto g_z_y_z_yy = cbuffer.data(pd_geom_11_off + 141 * ccomps * dcomps);

            auto g_z_y_z_yz = cbuffer.data(pd_geom_11_off + 142 * ccomps * dcomps);

            auto g_z_y_z_zz = cbuffer.data(pd_geom_11_off + 143 * ccomps * dcomps);

            auto g_z_z_x_xx = cbuffer.data(pd_geom_11_off + 144 * ccomps * dcomps);

            auto g_z_z_x_xy = cbuffer.data(pd_geom_11_off + 145 * ccomps * dcomps);

            auto g_z_z_x_xz = cbuffer.data(pd_geom_11_off + 146 * ccomps * dcomps);

            auto g_z_z_y_xx = cbuffer.data(pd_geom_11_off + 150 * ccomps * dcomps);

            auto g_z_z_y_xy = cbuffer.data(pd_geom_11_off + 151 * ccomps * dcomps);

            auto g_z_z_y_xz = cbuffer.data(pd_geom_11_off + 152 * ccomps * dcomps);

            auto g_z_z_y_yy = cbuffer.data(pd_geom_11_off + 153 * ccomps * dcomps);

            auto g_z_z_y_yz = cbuffer.data(pd_geom_11_off + 154 * ccomps * dcomps);

            auto g_z_z_z_xx = cbuffer.data(pd_geom_11_off + 156 * ccomps * dcomps);

            auto g_z_z_z_xy = cbuffer.data(pd_geom_11_off + 157 * ccomps * dcomps);

            auto g_z_z_z_xz = cbuffer.data(pd_geom_11_off + 158 * ccomps * dcomps);

            auto g_z_z_z_yy = cbuffer.data(pd_geom_11_off + 159 * ccomps * dcomps);

            auto g_z_z_z_yz = cbuffer.data(pd_geom_11_off + 160 * ccomps * dcomps);

            auto g_z_z_z_zz = cbuffer.data(pd_geom_11_off + 161 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dpxx

            const auto dp_geom_11_off = idx_geom_11_dpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_x_xx_x = cbuffer.data(dp_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_y = cbuffer.data(dp_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_z = cbuffer.data(dp_geom_11_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_x_0_x_x, g_x_0_x_y, g_x_0_x_z, g_x_x_x_x, g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_y, g_x_x_x_z, g_x_x_xx_x, g_x_x_xx_y, g_x_x_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xx_x[k] = -g_0_x_x_x[k] + g_x_0_x_x[k] - g_x_x_x_x[k] * ab_x + g_x_x_x_xx[k];

                g_x_x_xx_y[k] = -g_0_x_x_y[k] + g_x_0_x_y[k] - g_x_x_x_y[k] * ab_x + g_x_x_x_xy[k];

                g_x_x_xx_z[k] = -g_0_x_x_z[k] + g_x_0_x_z[k] - g_x_x_x_z[k] * ab_x + g_x_x_x_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_x_xy_x = cbuffer.data(dp_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xy_y = cbuffer.data(dp_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xy_z = cbuffer.data(dp_geom_11_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_xy, g_x_x_x_y, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_z, g_x_x_xy_x, g_x_x_xy_y, g_x_x_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xy_x[k] = -g_x_x_x_x[k] * ab_y + g_x_x_x_xy[k];

                g_x_x_xy_y[k] = -g_x_x_x_y[k] * ab_y + g_x_x_x_yy[k];

                g_x_x_xy_z[k] = -g_x_x_x_z[k] * ab_y + g_x_x_x_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_x_xz_x = cbuffer.data(dp_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xz_y = cbuffer.data(dp_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xz_z = cbuffer.data(dp_geom_11_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_xz, g_x_x_x_y, g_x_x_x_yz, g_x_x_x_z, g_x_x_x_zz, g_x_x_xz_x, g_x_x_xz_y, g_x_x_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xz_x[k] = -g_x_x_x_x[k] * ab_z + g_x_x_x_xz[k];

                g_x_x_xz_y[k] = -g_x_x_x_y[k] * ab_z + g_x_x_x_yz[k];

                g_x_x_xz_z[k] = -g_x_x_x_z[k] * ab_z + g_x_x_x_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_x_yy_x = cbuffer.data(dp_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_yy_y = cbuffer.data(dp_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_yy_z = cbuffer.data(dp_geom_11_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_xy, g_x_x_y_y, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_z, g_x_x_yy_x, g_x_x_yy_y, g_x_x_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yy_x[k] = -g_x_x_y_x[k] * ab_y + g_x_x_y_xy[k];

                g_x_x_yy_y[k] = -g_x_x_y_y[k] * ab_y + g_x_x_y_yy[k];

                g_x_x_yy_z[k] = -g_x_x_y_z[k] * ab_y + g_x_x_y_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_x_yz_x = cbuffer.data(dp_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_yz_y = cbuffer.data(dp_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_yz_z = cbuffer.data(dp_geom_11_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yz_x, g_x_x_yz_y, g_x_x_yz_z, g_x_x_z_x, g_x_x_z_xy, g_x_x_z_y, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yz_x[k] = -g_x_x_z_x[k] * ab_y + g_x_x_z_xy[k];

                g_x_x_yz_y[k] = -g_x_x_z_y[k] * ab_y + g_x_x_z_yy[k];

                g_x_x_yz_z[k] = -g_x_x_z_z[k] * ab_y + g_x_x_z_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_x_zz_x = cbuffer.data(dp_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_zz_y = cbuffer.data(dp_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_zz_z = cbuffer.data(dp_geom_11_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_xz, g_x_x_z_y, g_x_x_z_yz, g_x_x_z_z, g_x_x_z_zz, g_x_x_zz_x, g_x_x_zz_y, g_x_x_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zz_x[k] = -g_x_x_z_x[k] * ab_z + g_x_x_z_xz[k];

                g_x_x_zz_y[k] = -g_x_x_z_y[k] * ab_z + g_x_x_z_yz[k];

                g_x_x_zz_z[k] = -g_x_x_z_z[k] * ab_z + g_x_x_z_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_y_xx_x = cbuffer.data(dp_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_y_xx_y = cbuffer.data(dp_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_y_xx_z = cbuffer.data(dp_geom_11_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_x_y_x_x, g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_y, g_x_y_x_z, g_x_y_xx_x, g_x_y_xx_y, g_x_y_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xx_x[k] = -g_0_y_x_x[k] - g_x_y_x_x[k] * ab_x + g_x_y_x_xx[k];

                g_x_y_xx_y[k] = -g_0_y_x_y[k] - g_x_y_x_y[k] * ab_x + g_x_y_x_xy[k];

                g_x_y_xx_z[k] = -g_0_y_x_z[k] - g_x_y_x_z[k] * ab_x + g_x_y_x_xz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_y_xy_x = cbuffer.data(dp_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_y_xy_y = cbuffer.data(dp_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_y_xy_z = cbuffer.data(dp_geom_11_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_x_y_xy_x, g_x_y_xy_y, g_x_y_xy_z, g_x_y_y_x, g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_y, g_x_y_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xy_x[k] = -g_0_y_y_x[k] - g_x_y_y_x[k] * ab_x + g_x_y_y_xx[k];

                g_x_y_xy_y[k] = -g_0_y_y_y[k] - g_x_y_y_y[k] * ab_x + g_x_y_y_xy[k];

                g_x_y_xy_z[k] = -g_0_y_y_z[k] - g_x_y_y_z[k] * ab_x + g_x_y_y_xz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_y_xz_x = cbuffer.data(dp_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_y_xz_y = cbuffer.data(dp_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_y_xz_z = cbuffer.data(dp_geom_11_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_xz, g_x_y_x_y, g_x_y_x_yz, g_x_y_x_z, g_x_y_x_zz, g_x_y_xz_x, g_x_y_xz_y, g_x_y_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xz_x[k] = -g_x_y_x_x[k] * ab_z + g_x_y_x_xz[k];

                g_x_y_xz_y[k] = -g_x_y_x_y[k] * ab_z + g_x_y_x_yz[k];

                g_x_y_xz_z[k] = -g_x_y_x_z[k] * ab_z + g_x_y_x_zz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_x_y_yy_x = cbuffer.data(dp_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_y_yy_y = cbuffer.data(dp_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_y_yy_z = cbuffer.data(dp_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_x, g_x_0_y_y, g_x_0_y_z, g_x_y_y_x, g_x_y_y_xy, g_x_y_y_y, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_z, g_x_y_yy_x, g_x_y_yy_y, g_x_y_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yy_x[k] = g_x_0_y_x[k] - g_x_y_y_x[k] * ab_y + g_x_y_y_xy[k];

                g_x_y_yy_y[k] = g_x_0_y_y[k] - g_x_y_y_y[k] * ab_y + g_x_y_y_yy[k];

                g_x_y_yy_z[k] = g_x_0_y_z[k] - g_x_y_y_z[k] * ab_y + g_x_y_y_yz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_x_y_yz_x = cbuffer.data(dp_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_y_yz_y = cbuffer.data(dp_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_y_yz_z = cbuffer.data(dp_geom_11_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_xz, g_x_y_y_y, g_x_y_y_yz, g_x_y_y_z, g_x_y_y_zz, g_x_y_yz_x, g_x_y_yz_y, g_x_y_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yz_x[k] = -g_x_y_y_x[k] * ab_z + g_x_y_y_xz[k];

                g_x_y_yz_y[k] = -g_x_y_y_y[k] * ab_z + g_x_y_y_yz[k];

                g_x_y_yz_z[k] = -g_x_y_y_z[k] * ab_z + g_x_y_y_zz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_x_y_zz_x = cbuffer.data(dp_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_y_zz_y = cbuffer.data(dp_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_zz_z = cbuffer.data(dp_geom_11_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_xz, g_x_y_z_y, g_x_y_z_yz, g_x_y_z_z, g_x_y_z_zz, g_x_y_zz_x, g_x_y_zz_y, g_x_y_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zz_x[k] = -g_x_y_z_x[k] * ab_z + g_x_y_z_xz[k];

                g_x_y_zz_y[k] = -g_x_y_z_y[k] * ab_z + g_x_y_z_yz[k];

                g_x_y_zz_z[k] = -g_x_y_z_z[k] * ab_z + g_x_y_z_zz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_x_z_xx_x = cbuffer.data(dp_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_z_xx_y = cbuffer.data(dp_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_z_xx_z = cbuffer.data(dp_geom_11_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_x_z_x_x, g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_y, g_x_z_x_z, g_x_z_xx_x, g_x_z_xx_y, g_x_z_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xx_x[k] = -g_0_z_x_x[k] - g_x_z_x_x[k] * ab_x + g_x_z_x_xx[k];

                g_x_z_xx_y[k] = -g_0_z_x_y[k] - g_x_z_x_y[k] * ab_x + g_x_z_x_xy[k];

                g_x_z_xx_z[k] = -g_0_z_x_z[k] - g_x_z_x_z[k] * ab_x + g_x_z_x_xz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_x_z_xy_x = cbuffer.data(dp_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_z_xy_y = cbuffer.data(dp_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_z_xy_z = cbuffer.data(dp_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_x_x, g_x_z_x_xy, g_x_z_x_y, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_z, g_x_z_xy_x, g_x_z_xy_y, g_x_z_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xy_x[k] = -g_x_z_x_x[k] * ab_y + g_x_z_x_xy[k];

                g_x_z_xy_y[k] = -g_x_z_x_y[k] * ab_y + g_x_z_x_yy[k];

                g_x_z_xy_z[k] = -g_x_z_x_z[k] * ab_y + g_x_z_x_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_x_z_xz_x = cbuffer.data(dp_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_z_xz_y = cbuffer.data(dp_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_z_xz_z = cbuffer.data(dp_geom_11_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_x_z_xz_x, g_x_z_xz_y, g_x_z_xz_z, g_x_z_z_x, g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_y, g_x_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xz_x[k] = -g_0_z_z_x[k] - g_x_z_z_x[k] * ab_x + g_x_z_z_xx[k];

                g_x_z_xz_y[k] = -g_0_z_z_y[k] - g_x_z_z_y[k] * ab_x + g_x_z_z_xy[k];

                g_x_z_xz_z[k] = -g_0_z_z_z[k] - g_x_z_z_z[k] * ab_x + g_x_z_z_xz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_x_z_yy_x = cbuffer.data(dp_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_z_yy_y = cbuffer.data(dp_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_z_yy_z = cbuffer.data(dp_geom_11_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_y_x, g_x_z_y_xy, g_x_z_y_y, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_z, g_x_z_yy_x, g_x_z_yy_y, g_x_z_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yy_x[k] = -g_x_z_y_x[k] * ab_y + g_x_z_y_xy[k];

                g_x_z_yy_y[k] = -g_x_z_y_y[k] * ab_y + g_x_z_y_yy[k];

                g_x_z_yy_z[k] = -g_x_z_y_z[k] * ab_y + g_x_z_y_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_x_z_yz_x = cbuffer.data(dp_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_z_yz_y = cbuffer.data(dp_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_z_yz_z = cbuffer.data(dp_geom_11_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yz_x, g_x_z_yz_y, g_x_z_yz_z, g_x_z_z_x, g_x_z_z_xy, g_x_z_z_y, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yz_x[k] = -g_x_z_z_x[k] * ab_y + g_x_z_z_xy[k];

                g_x_z_yz_y[k] = -g_x_z_z_y[k] * ab_y + g_x_z_z_yy[k];

                g_x_z_yz_z[k] = -g_x_z_z_z[k] * ab_y + g_x_z_z_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_x_z_zz_x = cbuffer.data(dp_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_z_zz_y = cbuffer.data(dp_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_z_zz_z = cbuffer.data(dp_geom_11_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_x, g_x_0_z_y, g_x_0_z_z, g_x_z_z_x, g_x_z_z_xz, g_x_z_z_y, g_x_z_z_yz, g_x_z_z_z, g_x_z_z_zz, g_x_z_zz_x, g_x_z_zz_y, g_x_z_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zz_x[k] = g_x_0_z_x[k] - g_x_z_z_x[k] * ab_z + g_x_z_z_xz[k];

                g_x_z_zz_y[k] = g_x_0_z_y[k] - g_x_z_z_y[k] * ab_z + g_x_z_z_yz[k];

                g_x_z_zz_z[k] = g_x_0_z_z[k] - g_x_z_z_z[k] * ab_z + g_x_z_z_zz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_y_x_xx_x = cbuffer.data(dp_geom_11_off + 54 * ccomps * dcomps);

            auto g_y_x_xx_y = cbuffer.data(dp_geom_11_off + 55 * ccomps * dcomps);

            auto g_y_x_xx_z = cbuffer.data(dp_geom_11_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_x, g_y_0_x_y, g_y_0_x_z, g_y_x_x_x, g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_y, g_y_x_x_z, g_y_x_xx_x, g_y_x_xx_y, g_y_x_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xx_x[k] = g_y_0_x_x[k] - g_y_x_x_x[k] * ab_x + g_y_x_x_xx[k];

                g_y_x_xx_y[k] = g_y_0_x_y[k] - g_y_x_x_y[k] * ab_x + g_y_x_x_xy[k];

                g_y_x_xx_z[k] = g_y_0_x_z[k] - g_y_x_x_z[k] * ab_x + g_y_x_x_xz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_y_x_xy_x = cbuffer.data(dp_geom_11_off + 57 * ccomps * dcomps);

            auto g_y_x_xy_y = cbuffer.data(dp_geom_11_off + 58 * ccomps * dcomps);

            auto g_y_x_xy_z = cbuffer.data(dp_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_x, g_y_0_y_y, g_y_0_y_z, g_y_x_xy_x, g_y_x_xy_y, g_y_x_xy_z, g_y_x_y_x, g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_y, g_y_x_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xy_x[k] = g_y_0_y_x[k] - g_y_x_y_x[k] * ab_x + g_y_x_y_xx[k];

                g_y_x_xy_y[k] = g_y_0_y_y[k] - g_y_x_y_y[k] * ab_x + g_y_x_y_xy[k];

                g_y_x_xy_z[k] = g_y_0_y_z[k] - g_y_x_y_z[k] * ab_x + g_y_x_y_xz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_y_x_xz_x = cbuffer.data(dp_geom_11_off + 60 * ccomps * dcomps);

            auto g_y_x_xz_y = cbuffer.data(dp_geom_11_off + 61 * ccomps * dcomps);

            auto g_y_x_xz_z = cbuffer.data(dp_geom_11_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_xz, g_y_x_x_y, g_y_x_x_yz, g_y_x_x_z, g_y_x_x_zz, g_y_x_xz_x, g_y_x_xz_y, g_y_x_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xz_x[k] = -g_y_x_x_x[k] * ab_z + g_y_x_x_xz[k];

                g_y_x_xz_y[k] = -g_y_x_x_y[k] * ab_z + g_y_x_x_yz[k];

                g_y_x_xz_z[k] = -g_y_x_x_z[k] * ab_z + g_y_x_x_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_y_x_yy_x = cbuffer.data(dp_geom_11_off + 63 * ccomps * dcomps);

            auto g_y_x_yy_y = cbuffer.data(dp_geom_11_off + 64 * ccomps * dcomps);

            auto g_y_x_yy_z = cbuffer.data(dp_geom_11_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_y_x_y_x, g_y_x_y_xy, g_y_x_y_y, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_z, g_y_x_yy_x, g_y_x_yy_y, g_y_x_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yy_x[k] = -g_0_x_y_x[k] - g_y_x_y_x[k] * ab_y + g_y_x_y_xy[k];

                g_y_x_yy_y[k] = -g_0_x_y_y[k] - g_y_x_y_y[k] * ab_y + g_y_x_y_yy[k];

                g_y_x_yy_z[k] = -g_0_x_y_z[k] - g_y_x_y_z[k] * ab_y + g_y_x_y_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_y_x_yz_x = cbuffer.data(dp_geom_11_off + 66 * ccomps * dcomps);

            auto g_y_x_yz_y = cbuffer.data(dp_geom_11_off + 67 * ccomps * dcomps);

            auto g_y_x_yz_z = cbuffer.data(dp_geom_11_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_xz, g_y_x_y_y, g_y_x_y_yz, g_y_x_y_z, g_y_x_y_zz, g_y_x_yz_x, g_y_x_yz_y, g_y_x_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yz_x[k] = -g_y_x_y_x[k] * ab_z + g_y_x_y_xz[k];

                g_y_x_yz_y[k] = -g_y_x_y_y[k] * ab_z + g_y_x_y_yz[k];

                g_y_x_yz_z[k] = -g_y_x_y_z[k] * ab_z + g_y_x_y_zz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_y_x_zz_x = cbuffer.data(dp_geom_11_off + 69 * ccomps * dcomps);

            auto g_y_x_zz_y = cbuffer.data(dp_geom_11_off + 70 * ccomps * dcomps);

            auto g_y_x_zz_z = cbuffer.data(dp_geom_11_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_xz, g_y_x_z_y, g_y_x_z_yz, g_y_x_z_z, g_y_x_z_zz, g_y_x_zz_x, g_y_x_zz_y, g_y_x_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zz_x[k] = -g_y_x_z_x[k] * ab_z + g_y_x_z_xz[k];

                g_y_x_zz_y[k] = -g_y_x_z_y[k] * ab_z + g_y_x_z_yz[k];

                g_y_x_zz_z[k] = -g_y_x_z_z[k] * ab_z + g_y_x_z_zz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_y_y_xx_x = cbuffer.data(dp_geom_11_off + 72 * ccomps * dcomps);

            auto g_y_y_xx_y = cbuffer.data(dp_geom_11_off + 73 * ccomps * dcomps);

            auto g_y_y_xx_z = cbuffer.data(dp_geom_11_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_x_x, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_y, g_y_y_x_z, g_y_y_xx_x, g_y_y_xx_y, g_y_y_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xx_x[k] = -g_y_y_x_x[k] * ab_x + g_y_y_x_xx[k];

                g_y_y_xx_y[k] = -g_y_y_x_y[k] * ab_x + g_y_y_x_xy[k];

                g_y_y_xx_z[k] = -g_y_y_x_z[k] * ab_x + g_y_y_x_xz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_y_y_xy_x = cbuffer.data(dp_geom_11_off + 75 * ccomps * dcomps);

            auto g_y_y_xy_y = cbuffer.data(dp_geom_11_off + 76 * ccomps * dcomps);

            auto g_y_y_xy_z = cbuffer.data(dp_geom_11_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xy_x, g_y_y_xy_y, g_y_y_xy_z, g_y_y_y_x, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_y, g_y_y_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xy_x[k] = -g_y_y_y_x[k] * ab_x + g_y_y_y_xx[k];

                g_y_y_xy_y[k] = -g_y_y_y_y[k] * ab_x + g_y_y_y_xy[k];

                g_y_y_xy_z[k] = -g_y_y_y_z[k] * ab_x + g_y_y_y_xz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_y_y_xz_x = cbuffer.data(dp_geom_11_off + 78 * ccomps * dcomps);

            auto g_y_y_xz_y = cbuffer.data(dp_geom_11_off + 79 * ccomps * dcomps);

            auto g_y_y_xz_z = cbuffer.data(dp_geom_11_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xz_x, g_y_y_xz_y, g_y_y_xz_z, g_y_y_z_x, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_y, g_y_y_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xz_x[k] = -g_y_y_z_x[k] * ab_x + g_y_y_z_xx[k];

                g_y_y_xz_y[k] = -g_y_y_z_y[k] * ab_x + g_y_y_z_xy[k];

                g_y_y_xz_z[k] = -g_y_y_z_z[k] * ab_x + g_y_y_z_xz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_y_y_yy_x = cbuffer.data(dp_geom_11_off + 81 * ccomps * dcomps);

            auto g_y_y_yy_y = cbuffer.data(dp_geom_11_off + 82 * ccomps * dcomps);

            auto g_y_y_yy_z = cbuffer.data(dp_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_y_0_y_x, g_y_0_y_y, g_y_0_y_z, g_y_y_y_x, g_y_y_y_xy, g_y_y_y_y, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_z, g_y_y_yy_x, g_y_y_yy_y, g_y_y_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yy_x[k] = -g_0_y_y_x[k] + g_y_0_y_x[k] - g_y_y_y_x[k] * ab_y + g_y_y_y_xy[k];

                g_y_y_yy_y[k] = -g_0_y_y_y[k] + g_y_0_y_y[k] - g_y_y_y_y[k] * ab_y + g_y_y_y_yy[k];

                g_y_y_yy_z[k] = -g_0_y_y_z[k] + g_y_0_y_z[k] - g_y_y_y_z[k] * ab_y + g_y_y_y_yz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_y_y_yz_x = cbuffer.data(dp_geom_11_off + 84 * ccomps * dcomps);

            auto g_y_y_yz_y = cbuffer.data(dp_geom_11_off + 85 * ccomps * dcomps);

            auto g_y_y_yz_z = cbuffer.data(dp_geom_11_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_y_x, g_y_y_y_xz, g_y_y_y_y, g_y_y_y_yz, g_y_y_y_z, g_y_y_y_zz, g_y_y_yz_x, g_y_y_yz_y, g_y_y_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yz_x[k] = -g_y_y_y_x[k] * ab_z + g_y_y_y_xz[k];

                g_y_y_yz_y[k] = -g_y_y_y_y[k] * ab_z + g_y_y_y_yz[k];

                g_y_y_yz_z[k] = -g_y_y_y_z[k] * ab_z + g_y_y_y_zz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_y_y_zz_x = cbuffer.data(dp_geom_11_off + 87 * ccomps * dcomps);

            auto g_y_y_zz_y = cbuffer.data(dp_geom_11_off + 88 * ccomps * dcomps);

            auto g_y_y_zz_z = cbuffer.data(dp_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_z_x, g_y_y_z_xz, g_y_y_z_y, g_y_y_z_yz, g_y_y_z_z, g_y_y_z_zz, g_y_y_zz_x, g_y_y_zz_y, g_y_y_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zz_x[k] = -g_y_y_z_x[k] * ab_z + g_y_y_z_xz[k];

                g_y_y_zz_y[k] = -g_y_y_z_y[k] * ab_z + g_y_y_z_yz[k];

                g_y_y_zz_z[k] = -g_y_y_z_z[k] * ab_z + g_y_y_z_zz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_y_z_xx_x = cbuffer.data(dp_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_z_xx_y = cbuffer.data(dp_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_z_xx_z = cbuffer.data(dp_geom_11_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_x_x, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_y, g_y_z_x_z, g_y_z_xx_x, g_y_z_xx_y, g_y_z_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xx_x[k] = -g_y_z_x_x[k] * ab_x + g_y_z_x_xx[k];

                g_y_z_xx_y[k] = -g_y_z_x_y[k] * ab_x + g_y_z_x_xy[k];

                g_y_z_xx_z[k] = -g_y_z_x_z[k] * ab_x + g_y_z_x_xz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_y_z_xy_x = cbuffer.data(dp_geom_11_off + 93 * ccomps * dcomps);

            auto g_y_z_xy_y = cbuffer.data(dp_geom_11_off + 94 * ccomps * dcomps);

            auto g_y_z_xy_z = cbuffer.data(dp_geom_11_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xy_x, g_y_z_xy_y, g_y_z_xy_z, g_y_z_y_x, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_y, g_y_z_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xy_x[k] = -g_y_z_y_x[k] * ab_x + g_y_z_y_xx[k];

                g_y_z_xy_y[k] = -g_y_z_y_y[k] * ab_x + g_y_z_y_xy[k];

                g_y_z_xy_z[k] = -g_y_z_y_z[k] * ab_x + g_y_z_y_xz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_y_z_xz_x = cbuffer.data(dp_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_z_xz_y = cbuffer.data(dp_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_z_xz_z = cbuffer.data(dp_geom_11_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xz_x, g_y_z_xz_y, g_y_z_xz_z, g_y_z_z_x, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_y, g_y_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xz_x[k] = -g_y_z_z_x[k] * ab_x + g_y_z_z_xx[k];

                g_y_z_xz_y[k] = -g_y_z_z_y[k] * ab_x + g_y_z_z_xy[k];

                g_y_z_xz_z[k] = -g_y_z_z_z[k] * ab_x + g_y_z_z_xz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_y_z_yy_x = cbuffer.data(dp_geom_11_off + 99 * ccomps * dcomps);

            auto g_y_z_yy_y = cbuffer.data(dp_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_z_yy_z = cbuffer.data(dp_geom_11_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_y_z_y_x, g_y_z_y_xy, g_y_z_y_y, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_z, g_y_z_yy_x, g_y_z_yy_y, g_y_z_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yy_x[k] = -g_0_z_y_x[k] - g_y_z_y_x[k] * ab_y + g_y_z_y_xy[k];

                g_y_z_yy_y[k] = -g_0_z_y_y[k] - g_y_z_y_y[k] * ab_y + g_y_z_y_yy[k];

                g_y_z_yy_z[k] = -g_0_z_y_z[k] - g_y_z_y_z[k] * ab_y + g_y_z_y_yz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_y_z_yz_x = cbuffer.data(dp_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_z_yz_y = cbuffer.data(dp_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_z_yz_z = cbuffer.data(dp_geom_11_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_y_z_yz_x, g_y_z_yz_y, g_y_z_yz_z, g_y_z_z_x, g_y_z_z_xy, g_y_z_z_y, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yz_x[k] = -g_0_z_z_x[k] - g_y_z_z_x[k] * ab_y + g_y_z_z_xy[k];

                g_y_z_yz_y[k] = -g_0_z_z_y[k] - g_y_z_z_y[k] * ab_y + g_y_z_z_yy[k];

                g_y_z_yz_z[k] = -g_0_z_z_z[k] - g_y_z_z_z[k] * ab_y + g_y_z_z_yz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_y_z_zz_x = cbuffer.data(dp_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_z_zz_y = cbuffer.data(dp_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_z_zz_z = cbuffer.data(dp_geom_11_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_x, g_y_0_z_y, g_y_0_z_z, g_y_z_z_x, g_y_z_z_xz, g_y_z_z_y, g_y_z_z_yz, g_y_z_z_z, g_y_z_z_zz, g_y_z_zz_x, g_y_z_zz_y, g_y_z_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zz_x[k] = g_y_0_z_x[k] - g_y_z_z_x[k] * ab_z + g_y_z_z_xz[k];

                g_y_z_zz_y[k] = g_y_0_z_y[k] - g_y_z_z_y[k] * ab_z + g_y_z_z_yz[k];

                g_y_z_zz_z[k] = g_y_0_z_z[k] - g_y_z_z_z[k] * ab_z + g_y_z_z_zz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_z_x_xx_x = cbuffer.data(dp_geom_11_off + 108 * ccomps * dcomps);

            auto g_z_x_xx_y = cbuffer.data(dp_geom_11_off + 109 * ccomps * dcomps);

            auto g_z_x_xx_z = cbuffer.data(dp_geom_11_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_x, g_z_0_x_y, g_z_0_x_z, g_z_x_x_x, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_y, g_z_x_x_z, g_z_x_xx_x, g_z_x_xx_y, g_z_x_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xx_x[k] = g_z_0_x_x[k] - g_z_x_x_x[k] * ab_x + g_z_x_x_xx[k];

                g_z_x_xx_y[k] = g_z_0_x_y[k] - g_z_x_x_y[k] * ab_x + g_z_x_x_xy[k];

                g_z_x_xx_z[k] = g_z_0_x_z[k] - g_z_x_x_z[k] * ab_x + g_z_x_x_xz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_z_x_xy_x = cbuffer.data(dp_geom_11_off + 111 * ccomps * dcomps);

            auto g_z_x_xy_y = cbuffer.data(dp_geom_11_off + 112 * ccomps * dcomps);

            auto g_z_x_xy_z = cbuffer.data(dp_geom_11_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_x_x, g_z_x_x_xy, g_z_x_x_y, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_z, g_z_x_xy_x, g_z_x_xy_y, g_z_x_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xy_x[k] = -g_z_x_x_x[k] * ab_y + g_z_x_x_xy[k];

                g_z_x_xy_y[k] = -g_z_x_x_y[k] * ab_y + g_z_x_x_yy[k];

                g_z_x_xy_z[k] = -g_z_x_x_z[k] * ab_y + g_z_x_x_yz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_z_x_xz_x = cbuffer.data(dp_geom_11_off + 114 * ccomps * dcomps);

            auto g_z_x_xz_y = cbuffer.data(dp_geom_11_off + 115 * ccomps * dcomps);

            auto g_z_x_xz_z = cbuffer.data(dp_geom_11_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_x, g_z_0_z_y, g_z_0_z_z, g_z_x_xz_x, g_z_x_xz_y, g_z_x_xz_z, g_z_x_z_x, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_y, g_z_x_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xz_x[k] = g_z_0_z_x[k] - g_z_x_z_x[k] * ab_x + g_z_x_z_xx[k];

                g_z_x_xz_y[k] = g_z_0_z_y[k] - g_z_x_z_y[k] * ab_x + g_z_x_z_xy[k];

                g_z_x_xz_z[k] = g_z_0_z_z[k] - g_z_x_z_z[k] * ab_x + g_z_x_z_xz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_z_x_yy_x = cbuffer.data(dp_geom_11_off + 117 * ccomps * dcomps);

            auto g_z_x_yy_y = cbuffer.data(dp_geom_11_off + 118 * ccomps * dcomps);

            auto g_z_x_yy_z = cbuffer.data(dp_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_y_x, g_z_x_y_xy, g_z_x_y_y, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_z, g_z_x_yy_x, g_z_x_yy_y, g_z_x_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yy_x[k] = -g_z_x_y_x[k] * ab_y + g_z_x_y_xy[k];

                g_z_x_yy_y[k] = -g_z_x_y_y[k] * ab_y + g_z_x_y_yy[k];

                g_z_x_yy_z[k] = -g_z_x_y_z[k] * ab_y + g_z_x_y_yz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_z_x_yz_x = cbuffer.data(dp_geom_11_off + 120 * ccomps * dcomps);

            auto g_z_x_yz_y = cbuffer.data(dp_geom_11_off + 121 * ccomps * dcomps);

            auto g_z_x_yz_z = cbuffer.data(dp_geom_11_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yz_x, g_z_x_yz_y, g_z_x_yz_z, g_z_x_z_x, g_z_x_z_xy, g_z_x_z_y, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yz_x[k] = -g_z_x_z_x[k] * ab_y + g_z_x_z_xy[k];

                g_z_x_yz_y[k] = -g_z_x_z_y[k] * ab_y + g_z_x_z_yy[k];

                g_z_x_yz_z[k] = -g_z_x_z_z[k] * ab_y + g_z_x_z_yz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_z_x_zz_x = cbuffer.data(dp_geom_11_off + 123 * ccomps * dcomps);

            auto g_z_x_zz_y = cbuffer.data(dp_geom_11_off + 124 * ccomps * dcomps);

            auto g_z_x_zz_z = cbuffer.data(dp_geom_11_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_z_x_z_x, g_z_x_z_xz, g_z_x_z_y, g_z_x_z_yz, g_z_x_z_z, g_z_x_z_zz, g_z_x_zz_x, g_z_x_zz_y, g_z_x_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zz_x[k] = -g_0_x_z_x[k] - g_z_x_z_x[k] * ab_z + g_z_x_z_xz[k];

                g_z_x_zz_y[k] = -g_0_x_z_y[k] - g_z_x_z_y[k] * ab_z + g_z_x_z_yz[k];

                g_z_x_zz_z[k] = -g_0_x_z_z[k] - g_z_x_z_z[k] * ab_z + g_z_x_z_zz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_z_y_xx_x = cbuffer.data(dp_geom_11_off + 126 * ccomps * dcomps);

            auto g_z_y_xx_y = cbuffer.data(dp_geom_11_off + 127 * ccomps * dcomps);

            auto g_z_y_xx_z = cbuffer.data(dp_geom_11_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_x_x, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_y, g_z_y_x_z, g_z_y_xx_x, g_z_y_xx_y, g_z_y_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xx_x[k] = -g_z_y_x_x[k] * ab_x + g_z_y_x_xx[k];

                g_z_y_xx_y[k] = -g_z_y_x_y[k] * ab_x + g_z_y_x_xy[k];

                g_z_y_xx_z[k] = -g_z_y_x_z[k] * ab_x + g_z_y_x_xz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_z_y_xy_x = cbuffer.data(dp_geom_11_off + 129 * ccomps * dcomps);

            auto g_z_y_xy_y = cbuffer.data(dp_geom_11_off + 130 * ccomps * dcomps);

            auto g_z_y_xy_z = cbuffer.data(dp_geom_11_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xy_x, g_z_y_xy_y, g_z_y_xy_z, g_z_y_y_x, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_y, g_z_y_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xy_x[k] = -g_z_y_y_x[k] * ab_x + g_z_y_y_xx[k];

                g_z_y_xy_y[k] = -g_z_y_y_y[k] * ab_x + g_z_y_y_xy[k];

                g_z_y_xy_z[k] = -g_z_y_y_z[k] * ab_x + g_z_y_y_xz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_z_y_xz_x = cbuffer.data(dp_geom_11_off + 132 * ccomps * dcomps);

            auto g_z_y_xz_y = cbuffer.data(dp_geom_11_off + 133 * ccomps * dcomps);

            auto g_z_y_xz_z = cbuffer.data(dp_geom_11_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xz_x, g_z_y_xz_y, g_z_y_xz_z, g_z_y_z_x, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_y, g_z_y_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xz_x[k] = -g_z_y_z_x[k] * ab_x + g_z_y_z_xx[k];

                g_z_y_xz_y[k] = -g_z_y_z_y[k] * ab_x + g_z_y_z_xy[k];

                g_z_y_xz_z[k] = -g_z_y_z_z[k] * ab_x + g_z_y_z_xz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_z_y_yy_x = cbuffer.data(dp_geom_11_off + 135 * ccomps * dcomps);

            auto g_z_y_yy_y = cbuffer.data(dp_geom_11_off + 136 * ccomps * dcomps);

            auto g_z_y_yy_z = cbuffer.data(dp_geom_11_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_x, g_z_0_y_y, g_z_0_y_z, g_z_y_y_x, g_z_y_y_xy, g_z_y_y_y, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_z, g_z_y_yy_x, g_z_y_yy_y, g_z_y_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yy_x[k] = g_z_0_y_x[k] - g_z_y_y_x[k] * ab_y + g_z_y_y_xy[k];

                g_z_y_yy_y[k] = g_z_0_y_y[k] - g_z_y_y_y[k] * ab_y + g_z_y_y_yy[k];

                g_z_y_yy_z[k] = g_z_0_y_z[k] - g_z_y_y_z[k] * ab_y + g_z_y_y_yz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_z_y_yz_x = cbuffer.data(dp_geom_11_off + 138 * ccomps * dcomps);

            auto g_z_y_yz_y = cbuffer.data(dp_geom_11_off + 139 * ccomps * dcomps);

            auto g_z_y_yz_z = cbuffer.data(dp_geom_11_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_x, g_z_0_z_y, g_z_0_z_z, g_z_y_yz_x, g_z_y_yz_y, g_z_y_yz_z, g_z_y_z_x, g_z_y_z_xy, g_z_y_z_y, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yz_x[k] = g_z_0_z_x[k] - g_z_y_z_x[k] * ab_y + g_z_y_z_xy[k];

                g_z_y_yz_y[k] = g_z_0_z_y[k] - g_z_y_z_y[k] * ab_y + g_z_y_z_yy[k];

                g_z_y_yz_z[k] = g_z_0_z_z[k] - g_z_y_z_z[k] * ab_y + g_z_y_z_yz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_z_y_zz_x = cbuffer.data(dp_geom_11_off + 141 * ccomps * dcomps);

            auto g_z_y_zz_y = cbuffer.data(dp_geom_11_off + 142 * ccomps * dcomps);

            auto g_z_y_zz_z = cbuffer.data(dp_geom_11_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_z_y_z_x, g_z_y_z_xz, g_z_y_z_y, g_z_y_z_yz, g_z_y_z_z, g_z_y_z_zz, g_z_y_zz_x, g_z_y_zz_y, g_z_y_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zz_x[k] = -g_0_y_z_x[k] - g_z_y_z_x[k] * ab_z + g_z_y_z_xz[k];

                g_z_y_zz_y[k] = -g_0_y_z_y[k] - g_z_y_z_y[k] * ab_z + g_z_y_z_yz[k];

                g_z_y_zz_z[k] = -g_0_y_z_z[k] - g_z_y_z_z[k] * ab_z + g_z_y_z_zz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_z_z_xx_x = cbuffer.data(dp_geom_11_off + 144 * ccomps * dcomps);

            auto g_z_z_xx_y = cbuffer.data(dp_geom_11_off + 145 * ccomps * dcomps);

            auto g_z_z_xx_z = cbuffer.data(dp_geom_11_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_x_x, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_y, g_z_z_x_z, g_z_z_xx_x, g_z_z_xx_y, g_z_z_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xx_x[k] = -g_z_z_x_x[k] * ab_x + g_z_z_x_xx[k];

                g_z_z_xx_y[k] = -g_z_z_x_y[k] * ab_x + g_z_z_x_xy[k];

                g_z_z_xx_z[k] = -g_z_z_x_z[k] * ab_x + g_z_z_x_xz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_z_z_xy_x = cbuffer.data(dp_geom_11_off + 147 * ccomps * dcomps);

            auto g_z_z_xy_y = cbuffer.data(dp_geom_11_off + 148 * ccomps * dcomps);

            auto g_z_z_xy_z = cbuffer.data(dp_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xy_x, g_z_z_xy_y, g_z_z_xy_z, g_z_z_y_x, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_y, g_z_z_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xy_x[k] = -g_z_z_y_x[k] * ab_x + g_z_z_y_xx[k];

                g_z_z_xy_y[k] = -g_z_z_y_y[k] * ab_x + g_z_z_y_xy[k];

                g_z_z_xy_z[k] = -g_z_z_y_z[k] * ab_x + g_z_z_y_xz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_z_z_xz_x = cbuffer.data(dp_geom_11_off + 150 * ccomps * dcomps);

            auto g_z_z_xz_y = cbuffer.data(dp_geom_11_off + 151 * ccomps * dcomps);

            auto g_z_z_xz_z = cbuffer.data(dp_geom_11_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xz_x, g_z_z_xz_y, g_z_z_xz_z, g_z_z_z_x, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_y, g_z_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xz_x[k] = -g_z_z_z_x[k] * ab_x + g_z_z_z_xx[k];

                g_z_z_xz_y[k] = -g_z_z_z_y[k] * ab_x + g_z_z_z_xy[k];

                g_z_z_xz_z[k] = -g_z_z_z_z[k] * ab_x + g_z_z_z_xz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_z_z_yy_x = cbuffer.data(dp_geom_11_off + 153 * ccomps * dcomps);

            auto g_z_z_yy_y = cbuffer.data(dp_geom_11_off + 154 * ccomps * dcomps);

            auto g_z_z_yy_z = cbuffer.data(dp_geom_11_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_y_x, g_z_z_y_xy, g_z_z_y_y, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_z, g_z_z_yy_x, g_z_z_yy_y, g_z_z_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yy_x[k] = -g_z_z_y_x[k] * ab_y + g_z_z_y_xy[k];

                g_z_z_yy_y[k] = -g_z_z_y_y[k] * ab_y + g_z_z_y_yy[k];

                g_z_z_yy_z[k] = -g_z_z_y_z[k] * ab_y + g_z_z_y_yz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_z_z_yz_x = cbuffer.data(dp_geom_11_off + 156 * ccomps * dcomps);

            auto g_z_z_yz_y = cbuffer.data(dp_geom_11_off + 157 * ccomps * dcomps);

            auto g_z_z_yz_z = cbuffer.data(dp_geom_11_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yz_x, g_z_z_yz_y, g_z_z_yz_z, g_z_z_z_x, g_z_z_z_xy, g_z_z_z_y, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yz_x[k] = -g_z_z_z_x[k] * ab_y + g_z_z_z_xy[k];

                g_z_z_yz_y[k] = -g_z_z_z_y[k] * ab_y + g_z_z_z_yy[k];

                g_z_z_yz_z[k] = -g_z_z_z_z[k] * ab_y + g_z_z_z_yz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_z_z_zz_x = cbuffer.data(dp_geom_11_off + 159 * ccomps * dcomps);

            auto g_z_z_zz_y = cbuffer.data(dp_geom_11_off + 160 * ccomps * dcomps);

            auto g_z_z_zz_z = cbuffer.data(dp_geom_11_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_z_0_z_x, g_z_0_z_y, g_z_0_z_z, g_z_z_z_x, g_z_z_z_xz, g_z_z_z_y, g_z_z_z_yz, g_z_z_z_z, g_z_z_z_zz, g_z_z_zz_x, g_z_z_zz_y, g_z_z_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zz_x[k] = -g_0_z_z_x[k] + g_z_0_z_x[k] - g_z_z_z_x[k] * ab_z + g_z_z_z_xz[k];

                g_z_z_zz_y[k] = -g_0_z_z_y[k] + g_z_0_z_y[k] - g_z_z_z_y[k] * ab_z + g_z_z_z_yz[k];

                g_z_z_zz_z[k] = -g_0_z_z_z[k] + g_z_0_z_z[k] - g_z_z_z_z[k] * ab_z + g_z_z_z_zz[k];
            }
        }
    }
}

} // erirec namespace

