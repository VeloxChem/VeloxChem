#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxsp(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxsp,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxsp,
                                            const size_t idx_xxsd,
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
            /// Set up components of auxilary buffer : SSSP

            const auto sp_off = idx_xxsp + (i * bcomps + j) * 3;

            auto g_0_x = pbuffer.data(sp_off + 0);

            auto g_0_y = pbuffer.data(sp_off + 1);

            auto g_0_z = pbuffer.data(sp_off + 2);

            /// Set up components of auxilary buffer : SSSD

            const auto sd_off = idx_xxsd + (i * bcomps + j) * 6;

            auto g_0_xx = pbuffer.data(sd_off + 0);

            auto g_0_xy = pbuffer.data(sd_off + 1);

            auto g_0_xz = pbuffer.data(sd_off + 2);

            auto g_0_yy = pbuffer.data(sd_off + 3);

            auto g_0_yz = pbuffer.data(sd_off + 4);

            auto g_0_zz = pbuffer.data(sd_off + 5);

            /// set up bra offset for contr_buffer_xxsp

            const auto sp_geom_10_off = idx_geom_10_xxsp + (i * bcomps + j) * 3;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_0_x = cbuffer.data(sp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_0_y = cbuffer.data(sp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_0_z = cbuffer.data(sp_geom_10_off + 0 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_0_x, g_0_xx, g_0_xy, g_0_xz, g_0_y, g_0_z, g_x_0_0_x, g_x_0_0_y, g_x_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_x[k] = -g_0_x[k] * cd_x[k] + g_0_xx[k];

                g_x_0_0_y[k] = -g_0_y[k] * cd_x[k] + g_0_xy[k];

                g_x_0_0_z[k] = -g_0_z[k] * cd_x[k] + g_0_xz[k];
            }
            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_y_0_0_x = cbuffer.data(sp_geom_10_off + 3 * acomps * bcomps + 0);

            auto g_y_0_0_y = cbuffer.data(sp_geom_10_off + 3 * acomps * bcomps + 1);

            auto g_y_0_0_z = cbuffer.data(sp_geom_10_off + 3 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_y, g_0_x, g_0_xy, g_0_y, g_0_yy, g_0_yz, g_0_z, g_y_0_0_x, g_y_0_0_y, g_y_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_x[k] = -g_0_x[k] * cd_y[k] + g_0_xy[k];

                g_y_0_0_y[k] = -g_0_y[k] * cd_y[k] + g_0_yy[k];

                g_y_0_0_z[k] = -g_0_z[k] * cd_y[k] + g_0_yz[k];
            }
            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_z_0_0_x = cbuffer.data(sp_geom_10_off + 6 * acomps * bcomps + 0);

            auto g_z_0_0_y = cbuffer.data(sp_geom_10_off + 6 * acomps * bcomps + 1);

            auto g_z_0_0_z = cbuffer.data(sp_geom_10_off + 6 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_z, g_0_x, g_0_xz, g_0_y, g_0_yz, g_0_z, g_0_zz, g_z_0_0_x, g_z_0_0_y, g_z_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_x[k] = -g_0_x[k] * cd_z[k] + g_0_xz[k];

                g_z_0_0_y[k] = -g_0_y[k] * cd_z[k] + g_0_yz[k];

                g_z_0_0_z[k] = -g_0_z[k] * cd_z[k] + g_0_zz[k];
            }
        }
    }
}

} // erirec namespace

