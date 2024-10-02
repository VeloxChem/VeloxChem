#include "ElectronRepulsionGeomContrRecPDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec {  // erirec namespace

auto comp_bra_geom_hrr_electron_repulsion_pdxx(CSimdArray<double>&   cbuffer,
                                               const size_t          idx_geom_pdxx,
                                               const size_t          idx_geom_sdxx,
                                               const size_t          idx_geom_sfxx,
                                               const size_t          idx_sdxx,
                                               const TPoint<double>& r_ab,
                                               const int             c_angmom,
                                               const int             d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{
        c_angmom,
    });

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{
        d_angmom,
    });

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

            auto sd_geom_off = idx_geom_sdxx + i * dcomps + j;

            auto g_x_0_xx = cbuffer.data(sd_geom_off + 0 * ccomps * dcomps);

            auto g_x_0_xy = cbuffer.data(sd_geom_off + 1 * ccomps * dcomps);

            auto g_x_0_xz = cbuffer.data(sd_geom_off + 2 * ccomps * dcomps);

            auto g_x_0_yy = cbuffer.data(sd_geom_off + 3 * ccomps * dcomps);

            auto g_x_0_yz = cbuffer.data(sd_geom_off + 4 * ccomps * dcomps);

            auto g_x_0_zz = cbuffer.data(sd_geom_off + 5 * ccomps * dcomps);
            
            sd_geom_off += 6 * ccomps * dcomps;
            
            auto g_y_0_xx = cbuffer.data(sd_geom_off + 0 * ccomps * dcomps);

            auto g_y_0_xy = cbuffer.data(sd_geom_off + 1 * ccomps * dcomps);

            auto g_y_0_xz = cbuffer.data(sd_geom_off + 2 * ccomps * dcomps);

            auto g_y_0_yy = cbuffer.data(sd_geom_off + 3 * ccomps * dcomps);

            auto g_y_0_yz = cbuffer.data(sd_geom_off + 4 * ccomps * dcomps);

            auto g_y_0_zz = cbuffer.data(sd_geom_off + 5 * ccomps * dcomps);
            
            sd_geom_off += 6 * ccomps * dcomps;
            
            auto g_z_0_xx = cbuffer.data(sd_geom_off + 0 * ccomps * dcomps);

            auto g_z_0_xy = cbuffer.data(sd_geom_off + 1 * ccomps * dcomps);

            auto g_z_0_xz = cbuffer.data(sd_geom_off + 2 * ccomps * dcomps);

            auto g_z_0_yy = cbuffer.data(sd_geom_off + 3 * ccomps * dcomps);

            auto g_z_0_yz = cbuffer.data(sd_geom_off + 4 * ccomps * dcomps);

            auto g_z_0_zz = cbuffer.data(sd_geom_off + 5 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SFSS

            auto sf_geom_off = idx_geom_sfxx + i * dcomps + j;

            auto g_x_0_xxx = cbuffer.data(sf_geom_off + 0 * ccomps * dcomps);

            auto g_x_0_xxy = cbuffer.data(sf_geom_off + 1 * ccomps * dcomps);

            auto g_x_0_xxz = cbuffer.data(sf_geom_off + 2 * ccomps * dcomps);

            auto g_x_0_xyy = cbuffer.data(sf_geom_off + 3 * ccomps * dcomps);

            auto g_x_0_xyz = cbuffer.data(sf_geom_off + 4 * ccomps * dcomps);

            auto g_x_0_xzz = cbuffer.data(sf_geom_off + 5 * ccomps * dcomps);

            auto g_x_0_yyy = cbuffer.data(sf_geom_off + 6 * ccomps * dcomps);

            auto g_x_0_yyz = cbuffer.data(sf_geom_off + 7 * ccomps * dcomps);

            auto g_x_0_yzz = cbuffer.data(sf_geom_off + 8 * ccomps * dcomps);

            auto g_x_0_zzz = cbuffer.data(sf_geom_off + 9 * ccomps * dcomps);
            
            sf_geom_off += 10 * ccomps * dcomps;
            
            auto g_y_0_xxx = cbuffer.data(sf_geom_off + 0 * ccomps * dcomps);

            auto g_y_0_xxy = cbuffer.data(sf_geom_off + 1 * ccomps * dcomps);

            auto g_y_0_xxz = cbuffer.data(sf_geom_off + 2 * ccomps * dcomps);

            auto g_y_0_xyy = cbuffer.data(sf_geom_off + 3 * ccomps * dcomps);

            auto g_y_0_xyz = cbuffer.data(sf_geom_off + 4 * ccomps * dcomps);

            auto g_y_0_xzz = cbuffer.data(sf_geom_off + 5 * ccomps * dcomps);

            auto g_y_0_yyy = cbuffer.data(sf_geom_off + 6 * ccomps * dcomps);

            auto g_y_0_yyz = cbuffer.data(sf_geom_off + 7 * ccomps * dcomps);

            auto g_y_0_yzz = cbuffer.data(sf_geom_off + 8 * ccomps * dcomps);

            auto g_y_0_zzz = cbuffer.data(sf_geom_off + 9 * ccomps * dcomps);
            
            sf_geom_off += 10 * ccomps * dcomps;
            
            auto g_z_0_xxx = cbuffer.data(sf_geom_off + 0 * ccomps * dcomps);

            auto g_z_0_xxy = cbuffer.data(sf_geom_off + 1 * ccomps * dcomps);

            auto g_z_0_xxz = cbuffer.data(sf_geom_off + 2 * ccomps * dcomps);

            auto g_z_0_xyy = cbuffer.data(sf_geom_off + 3 * ccomps * dcomps);

            auto g_z_0_xyz = cbuffer.data(sf_geom_off + 4 * ccomps * dcomps);

            auto g_z_0_xzz = cbuffer.data(sf_geom_off + 5 * ccomps * dcomps);

            auto g_z_0_yyy = cbuffer.data(sf_geom_off + 6 * ccomps * dcomps);

            auto g_z_0_yyz = cbuffer.data(sf_geom_off + 7 * ccomps * dcomps);

            auto g_z_0_yzz = cbuffer.data(sf_geom_off + 8 * ccomps * dcomps);

            auto g_z_0_zzz = cbuffer.data(sf_geom_off + 9 * ccomps * dcomps);
            
            /// Set up components of auxilary buffer : SDSS

            auto sd_off = idx_sdxx + i * dcomps + j;

            auto g_0_xx = cbuffer.data(sd_off + 0 * ccomps * dcomps);

            auto g_0_xy = cbuffer.data(sd_off + 1 * ccomps * dcomps);

            auto g_0_xz = cbuffer.data(sd_off + 2 * ccomps * dcomps);

            auto g_0_yy = cbuffer.data(sd_off + 3 * ccomps * dcomps);

            auto g_0_yz = cbuffer.data(sd_off + 4 * ccomps * dcomps);

            auto g_0_zz = cbuffer.data(sd_off + 5 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pdxx

            auto pd_geom_off = idx_geom_pdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_x_xx = cbuffer.data(pd_geom_off + 0 * ccomps * dcomps);

            auto g_x_x_xy = cbuffer.data(pd_geom_off + 1 * ccomps * dcomps);

            auto g_x_x_xz = cbuffer.data(pd_geom_off + 2 * ccomps * dcomps);

            auto g_x_x_yy = cbuffer.data(pd_geom_off + 3 * ccomps * dcomps);

            auto g_x_x_yz = cbuffer.data(pd_geom_off + 4 * ccomps * dcomps);

            auto g_x_x_zz = cbuffer.data(pd_geom_off + 5 * ccomps * dcomps);

#pragma omp simd aligned(g_x_0_xx,      \
                             g_x_0_xxx, \
                             g_x_0_xxy, \
                             g_x_0_xxz, \
                             g_x_0_xy,  \
                             g_x_0_xyy, \
                             g_x_0_xyz, \
                             g_x_0_xz,  \
                             g_x_0_xzz, \
                             g_x_0_yy,  \
                             g_x_0_yz,  \
                             g_x_0_zz,  \
                             g_0_xx,    \
                             g_0_xy,    \
                             g_0_xz,    \
                             g_0_yy,    \
                             g_0_yz,    \
                             g_0_zz,    \
                             g_x_x_xx,  \
                             g_x_x_xy,  \
                             g_x_x_xz,  \
                             g_x_x_yy,  \
                             g_x_x_yz,  \
                             g_x_x_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xx[k] = -g_x_0_xx[k] * ab_x + g_x_0_xxx[k] - g_0_xx[k];

                g_x_x_xy[k] = -g_x_0_xy[k] * ab_x + g_x_0_xxy[k] - g_0_xy[k];

                g_x_x_xz[k] = -g_x_0_xz[k] * ab_x + g_x_0_xxz[k] - g_0_xz[k];

                g_x_x_yy[k] = -g_x_0_yy[k] * ab_x + g_x_0_xyy[k] - g_0_yy[k];

                g_x_x_yz[k] = -g_x_0_yz[k] * ab_x + g_x_0_xyz[k] - g_0_yz[k];

                g_x_x_zz[k] = -g_x_0_zz[k] * ab_x + g_x_0_xzz[k] - g_0_zz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_x_y_xx = cbuffer.data(pd_geom_off + 6 * ccomps * dcomps);

            auto g_x_y_xy = cbuffer.data(pd_geom_off + 7 * ccomps * dcomps);

            auto g_x_y_xz = cbuffer.data(pd_geom_off + 8 * ccomps * dcomps);

            auto g_x_y_yy = cbuffer.data(pd_geom_off + 9 * ccomps * dcomps);

            auto g_x_y_yz = cbuffer.data(pd_geom_off + 10 * ccomps * dcomps);

            auto g_x_y_zz = cbuffer.data(pd_geom_off + 11 * ccomps * dcomps);

#pragma omp simd aligned(g_x_0_xx,      \
                             g_x_0_xxy, \
                             g_x_0_xy,  \
                             g_x_0_xyy, \
                             g_x_0_xyz, \
                             g_x_0_xz,  \
                             g_x_0_yy,  \
                             g_x_0_yyy, \
                             g_x_0_yyz, \
                             g_x_0_yz,  \
                             g_x_0_yzz, \
                             g_x_0_zz,  \
                             g_x_y_xx,  \
                             g_x_y_xy,  \
                             g_x_y_xz,  \
                             g_x_y_yy,  \
                             g_x_y_yz,  \
                             g_x_y_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xx[k] = -g_x_0_xx[k] * ab_y + g_x_0_xxy[k];

                g_x_y_xy[k] = -g_x_0_xy[k] * ab_y + g_x_0_xyy[k];

                g_x_y_xz[k] = -g_x_0_xz[k] * ab_y + g_x_0_xyz[k];

                g_x_y_yy[k] = -g_x_0_yy[k] * ab_y + g_x_0_yyy[k];

                g_x_y_yz[k] = -g_x_0_yz[k] * ab_y + g_x_0_yyz[k];

                g_x_y_zz[k] = -g_x_0_zz[k] * ab_y + g_x_0_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_x_z_xx = cbuffer.data(pd_geom_off + 12 * ccomps * dcomps);

            auto g_x_z_xy = cbuffer.data(pd_geom_off + 13 * ccomps * dcomps);

            auto g_x_z_xz = cbuffer.data(pd_geom_off + 14 * ccomps * dcomps);

            auto g_x_z_yy = cbuffer.data(pd_geom_off + 15 * ccomps * dcomps);

            auto g_x_z_yz = cbuffer.data(pd_geom_off + 16 * ccomps * dcomps);

            auto g_x_z_zz = cbuffer.data(pd_geom_off + 17 * ccomps * dcomps);

#pragma omp simd aligned(g_x_0_xx,      \
                             g_x_0_xxz, \
                             g_x_0_xy,  \
                             g_x_0_xyz, \
                             g_x_0_xz,  \
                             g_x_0_xzz, \
                             g_x_0_yy,  \
                             g_x_0_yyz, \
                             g_x_0_yz,  \
                             g_x_0_yzz, \
                             g_x_0_zz,  \
                             g_x_0_zzz, \
                             g_x_z_xx,  \
                             g_x_z_xy,  \
                             g_x_z_xz,  \
                             g_x_z_yy,  \
                             g_x_z_yz,  \
                             g_x_z_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xx[k] = -g_x_0_xx[k] * ab_z + g_x_0_xxz[k];

                g_x_z_xy[k] = -g_x_0_xy[k] * ab_z + g_x_0_xyz[k];

                g_x_z_xz[k] = -g_x_0_xz[k] * ab_z + g_x_0_xzz[k];

                g_x_z_yy[k] = -g_x_0_yy[k] * ab_z + g_x_0_yyz[k];

                g_x_z_yz[k] = -g_x_0_yz[k] * ab_z + g_x_0_yzz[k];

                g_x_z_zz[k] = -g_x_0_zz[k] * ab_z + g_x_0_zzz[k];
            }
            
            // updated offset
            
            pd_geom_off += 18 * ccomps * dcomps;
            
            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_y_x_xx = cbuffer.data(pd_geom_off + 0 * ccomps * dcomps);

            auto g_y_x_xy = cbuffer.data(pd_geom_off + 1 * ccomps * dcomps);

            auto g_y_x_xz = cbuffer.data(pd_geom_off + 2 * ccomps * dcomps);

            auto g_y_x_yy = cbuffer.data(pd_geom_off + 3 * ccomps * dcomps);

            auto g_y_x_yz = cbuffer.data(pd_geom_off + 4 * ccomps * dcomps);

            auto g_y_x_zz = cbuffer.data(pd_geom_off + 5 * ccomps * dcomps);

#pragma omp simd aligned(g_y_0_xx,      \
                             g_y_0_xxx, \
                             g_y_0_xxy, \
                             g_y_0_xxz, \
                             g_y_0_xy,  \
                             g_y_0_xyy, \
                             g_y_0_xyz, \
                             g_y_0_xz,  \
                             g_y_0_xzz, \
                             g_y_0_yy,  \
                             g_y_0_yz,  \
                             g_y_0_zz,  \
                             g_y_x_xx,  \
                             g_y_x_xy,  \
                             g_y_x_xz,  \
                             g_y_x_yy,  \
                             g_y_x_yz,  \
                             g_y_x_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xx[k] = -g_y_0_xx[k] * ab_x + g_y_0_xxx[k];

                g_y_x_xy[k] = -g_y_0_xy[k] * ab_x + g_y_0_xxy[k];

                g_y_x_xz[k] = -g_y_0_xz[k] * ab_x + g_y_0_xxz[k];

                g_y_x_yy[k] = -g_y_0_yy[k] * ab_x + g_y_0_xyy[k];

                g_y_x_yz[k] = -g_y_0_yz[k] * ab_x + g_y_0_xyz[k];

                g_y_x_zz[k] = -g_y_0_zz[k] * ab_x + g_y_0_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_y_y_xx = cbuffer.data(pd_geom_off + 6 * ccomps * dcomps);

            auto g_y_y_xy = cbuffer.data(pd_geom_off + 7 * ccomps * dcomps);

            auto g_y_y_xz = cbuffer.data(pd_geom_off + 8 * ccomps * dcomps);

            auto g_y_y_yy = cbuffer.data(pd_geom_off + 9 * ccomps * dcomps);

            auto g_y_y_yz = cbuffer.data(pd_geom_off + 10 * ccomps * dcomps);

            auto g_y_y_zz = cbuffer.data(pd_geom_off + 11 * ccomps * dcomps);

#pragma omp simd aligned(g_y_0_xx,      \
                             g_y_0_xxy, \
                             g_y_0_xy,  \
                             g_y_0_xyy, \
                             g_y_0_xyz, \
                             g_y_0_xz,  \
                             g_y_0_yy,  \
                             g_y_0_yyy, \
                             g_y_0_yyz, \
                             g_y_0_yz,  \
                             g_y_0_yzz, \
                             g_y_0_zz,  \
                             g_0_xx,    \
                             g_0_xy,    \
                             g_0_xz,    \
                             g_0_yy,    \
                             g_0_yz,    \
                             g_0_zz,    \
                             g_y_y_xx,  \
                             g_y_y_xy,  \
                             g_y_y_xz,  \
                             g_y_y_yy,  \
                             g_y_y_yz,  \
                             g_y_y_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xx[k] = -g_y_0_xx[k] * ab_y + g_y_0_xxy[k] - g_0_xx[k];

                g_y_y_xy[k] = -g_y_0_xy[k] * ab_y + g_y_0_xyy[k] - g_0_xy[k];

                g_y_y_xz[k] = -g_y_0_xz[k] * ab_y + g_y_0_xyz[k] - g_0_xz[k];

                g_y_y_yy[k] = -g_y_0_yy[k] * ab_y + g_y_0_yyy[k] - g_0_yy[k];

                g_y_y_yz[k] = -g_y_0_yz[k] * ab_y + g_y_0_yyz[k] - g_0_yz[k];

                g_y_y_zz[k] = -g_y_0_zz[k] * ab_y + g_y_0_yzz[k] - g_0_zz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_y_z_xx = cbuffer.data(pd_geom_off + 12 * ccomps * dcomps);

            auto g_y_z_xy = cbuffer.data(pd_geom_off + 13 * ccomps * dcomps);

            auto g_y_z_xz = cbuffer.data(pd_geom_off + 14 * ccomps * dcomps);

            auto g_y_z_yy = cbuffer.data(pd_geom_off + 15 * ccomps * dcomps);

            auto g_y_z_yz = cbuffer.data(pd_geom_off + 16 * ccomps * dcomps);

            auto g_y_z_zz = cbuffer.data(pd_geom_off + 17 * ccomps * dcomps);

#pragma omp simd aligned(g_y_0_xx,      \
                             g_y_0_xxz, \
                             g_y_0_xy,  \
                             g_y_0_xyz, \
                             g_y_0_xz,  \
                             g_y_0_xzz, \
                             g_y_0_yy,  \
                             g_y_0_yyz, \
                             g_y_0_yz,  \
                             g_y_0_yzz, \
                             g_y_0_zz,  \
                             g_y_0_zzz, \
                             g_y_z_xx,  \
                             g_y_z_xy,  \
                             g_y_z_xz,  \
                             g_y_z_yy,  \
                             g_y_z_yz,  \
                             g_y_z_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xx[k] = -g_y_0_xx[k] * ab_z + g_y_0_xxz[k];

                g_y_z_xy[k] = -g_y_0_xy[k] * ab_z + g_y_0_xyz[k];

                g_y_z_xz[k] = -g_y_0_xz[k] * ab_z + g_y_0_xzz[k];

                g_y_z_yy[k] = -g_y_0_yy[k] * ab_z + g_y_0_yyz[k];

                g_y_z_yz[k] = -g_y_0_yz[k] * ab_z + g_y_0_yzz[k];

                g_y_z_zz[k] = -g_y_0_zz[k] * ab_z + g_y_0_zzz[k];
            }
            
            // updated offset
            
            pd_geom_off += 18 * ccomps * dcomps;
            
            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_z_x_xx = cbuffer.data(pd_geom_off + 0 * ccomps * dcomps);

            auto g_z_x_xy = cbuffer.data(pd_geom_off + 1 * ccomps * dcomps);

            auto g_z_x_xz = cbuffer.data(pd_geom_off + 2 * ccomps * dcomps);

            auto g_z_x_yy = cbuffer.data(pd_geom_off + 3 * ccomps * dcomps);

            auto g_z_x_yz = cbuffer.data(pd_geom_off + 4 * ccomps * dcomps);

            auto g_z_x_zz = cbuffer.data(pd_geom_off + 5 * ccomps * dcomps);

#pragma omp simd aligned(g_z_0_xx,      \
                             g_z_0_xxx, \
                             g_z_0_xxy, \
                             g_z_0_xxz, \
                             g_z_0_xy,  \
                             g_z_0_xyy, \
                             g_z_0_xyz, \
                             g_z_0_xz,  \
                             g_z_0_xzz, \
                             g_z_0_yy,  \
                             g_z_0_yz,  \
                             g_z_0_zz,  \
                             g_z_x_xx,  \
                             g_z_x_xy,  \
                             g_z_x_xz,  \
                             g_z_x_yy,  \
                             g_z_x_yz,  \
                             g_z_x_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xx[k] = -g_z_0_xx[k] * ab_x + g_z_0_xxx[k];

                g_z_x_xy[k] = -g_z_0_xy[k] * ab_x + g_z_0_xxy[k];

                g_z_x_xz[k] = -g_z_0_xz[k] * ab_x + g_z_0_xxz[k];

                g_z_x_yy[k] = -g_z_0_yy[k] * ab_x + g_z_0_xyy[k];

                g_z_x_yz[k] = -g_z_0_yz[k] * ab_x + g_z_0_xyz[k];

                g_z_x_zz[k] = -g_z_0_zz[k] * ab_x + g_z_0_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_z_y_xx = cbuffer.data(pd_geom_off + 6 * ccomps * dcomps);

            auto g_z_y_xy = cbuffer.data(pd_geom_off + 7 * ccomps * dcomps);

            auto g_z_y_xz = cbuffer.data(pd_geom_off + 8 * ccomps * dcomps);

            auto g_z_y_yy = cbuffer.data(pd_geom_off + 9 * ccomps * dcomps);

            auto g_z_y_yz = cbuffer.data(pd_geom_off + 10 * ccomps * dcomps);

            auto g_z_y_zz = cbuffer.data(pd_geom_off + 11 * ccomps * dcomps);

#pragma omp simd aligned(g_z_0_xx,      \
                             g_z_0_xxy, \
                             g_z_0_xy,  \
                             g_z_0_xyy, \
                             g_z_0_xyz, \
                             g_z_0_xz,  \
                             g_z_0_yy,  \
                             g_z_0_yyy, \
                             g_z_0_yyz, \
                             g_z_0_yz,  \
                             g_z_0_yzz, \
                             g_z_0_zz,  \
                             g_z_y_xx,  \
                             g_z_y_xy,  \
                             g_z_y_xz,  \
                             g_z_y_yy,  \
                             g_z_y_yz,  \
                             g_z_y_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xx[k] = -g_z_0_xx[k] * ab_y + g_z_0_xxy[k];

                g_z_y_xy[k] = -g_z_0_xy[k] * ab_y + g_z_0_xyy[k];

                g_z_y_xz[k] = -g_z_0_xz[k] * ab_y + g_z_0_xyz[k];

                g_z_y_yy[k] = -g_z_0_yy[k] * ab_y + g_z_0_yyy[k];

                g_z_y_yz[k] = -g_z_0_yz[k] * ab_y + g_z_0_yyz[k];

                g_z_y_zz[k] = -g_z_0_zz[k] * ab_y + g_z_0_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_z_z_xx = cbuffer.data(pd_geom_off + 12 * ccomps * dcomps);

            auto g_z_z_xy = cbuffer.data(pd_geom_off + 13 * ccomps * dcomps);

            auto g_z_z_xz = cbuffer.data(pd_geom_off + 14 * ccomps * dcomps);

            auto g_z_z_yy = cbuffer.data(pd_geom_off + 15 * ccomps * dcomps);

            auto g_z_z_yz = cbuffer.data(pd_geom_off + 16 * ccomps * dcomps);

            auto g_z_z_zz = cbuffer.data(pd_geom_off + 17 * ccomps * dcomps);

#pragma omp simd aligned(g_z_0_xx,      \
                             g_z_0_xxz, \
                             g_z_0_xy,  \
                             g_z_0_xyz, \
                             g_z_0_xz,  \
                             g_z_0_xzz, \
                             g_z_0_yy,  \
                             g_z_0_yyz, \
                             g_z_0_yz,  \
                             g_z_0_yzz, \
                             g_z_0_zz,  \
                             g_z_0_zzz, \
                             g_0_xx,    \
                             g_0_xy,    \
                             g_0_xz,    \
                             g_0_yy,    \
                             g_0_yz,    \
                             g_0_zz,    \
                             g_z_z_xx,  \
                             g_z_z_xy,  \
                             g_z_z_xz,  \
                             g_z_z_yy,  \
                             g_z_z_yz,  \
                             g_z_z_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xx[k] = -g_z_0_xx[k] * ab_z + g_z_0_xxz[k] - g_0_xx[k];

                g_z_z_xy[k] = -g_z_0_xy[k] * ab_z + g_z_0_xyz[k] - g_0_xy[k];

                g_z_z_xz[k] = -g_z_0_xz[k] * ab_z + g_z_0_xzz[k] - g_0_xz[k];

                g_z_z_yy[k] = -g_z_0_yy[k] * ab_z + g_z_0_yyz[k] - g_0_yy[k];

                g_z_z_yz[k] = -g_z_0_yz[k] * ab_z + g_z_0_yzz[k] - g_0_yz[k];

                g_z_z_zz[k] = -g_z_0_zz[k] * ab_z + g_z_0_zzz[k] - g_0_zz[k];
            }
        }
    }
}


}  // namespace erirec
