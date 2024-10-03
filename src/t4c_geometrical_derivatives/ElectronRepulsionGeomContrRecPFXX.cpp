#include "ElectronRepulsionGeomContrRecPFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec {  // erirec namespace

auto comp_bra_geom_hrr_electron_repulsion_pfxx(CSimdArray<double>&   cbuffer,
                                               const size_t          idx_geom_pfxx,
                                               const size_t          idx_geom_sfxx,
                                               const size_t          idx_geom_sgxx,
                                               const size_t          idx_sfxx,
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
            
            /// Set up components of auxilary buffer : SGSS

            auto sg_geom_off = idx_geom_sgxx + i * dcomps + j;

            auto g_x_0_xxxx = cbuffer.data(sg_geom_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxy = cbuffer.data(sg_geom_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxz = cbuffer.data(sg_geom_off + 2 * ccomps * dcomps);

            auto g_x_0_xxyy = cbuffer.data(sg_geom_off + 3 * ccomps * dcomps);

            auto g_x_0_xxyz = cbuffer.data(sg_geom_off + 4 * ccomps * dcomps);

            auto g_x_0_xxzz = cbuffer.data(sg_geom_off + 5 * ccomps * dcomps);

            auto g_x_0_xyyy = cbuffer.data(sg_geom_off + 6 * ccomps * dcomps);

            auto g_x_0_xyyz = cbuffer.data(sg_geom_off + 7 * ccomps * dcomps);

            auto g_x_0_xyzz = cbuffer.data(sg_geom_off + 8 * ccomps * dcomps);

            auto g_x_0_xzzz = cbuffer.data(sg_geom_off + 9 * ccomps * dcomps);

            auto g_x_0_yyyy = cbuffer.data(sg_geom_off + 10 * ccomps * dcomps);

            auto g_x_0_yyyz = cbuffer.data(sg_geom_off + 11 * ccomps * dcomps);

            auto g_x_0_yyzz = cbuffer.data(sg_geom_off + 12 * ccomps * dcomps);

            auto g_x_0_yzzz = cbuffer.data(sg_geom_off + 13 * ccomps * dcomps);

            auto g_x_0_zzzz = cbuffer.data(sg_geom_off + 14 * ccomps * dcomps);
            
            sg_geom_off += 15 * ccomps * dcomps;
            
            auto g_y_0_xxxx = cbuffer.data(sg_geom_off + 0 * ccomps * dcomps);

            auto g_y_0_xxxy = cbuffer.data(sg_geom_off + 1 * ccomps * dcomps);

            auto g_y_0_xxxz = cbuffer.data(sg_geom_off + 2 * ccomps * dcomps);

            auto g_y_0_xxyy = cbuffer.data(sg_geom_off + 3 * ccomps * dcomps);

            auto g_y_0_xxyz = cbuffer.data(sg_geom_off + 4 * ccomps * dcomps);

            auto g_y_0_xxzz = cbuffer.data(sg_geom_off + 5 * ccomps * dcomps);

            auto g_y_0_xyyy = cbuffer.data(sg_geom_off + 6 * ccomps * dcomps);

            auto g_y_0_xyyz = cbuffer.data(sg_geom_off + 7 * ccomps * dcomps);

            auto g_y_0_xyzz = cbuffer.data(sg_geom_off + 8 * ccomps * dcomps);

            auto g_y_0_xzzz = cbuffer.data(sg_geom_off + 9 * ccomps * dcomps);

            auto g_y_0_yyyy = cbuffer.data(sg_geom_off + 10 * ccomps * dcomps);

            auto g_y_0_yyyz = cbuffer.data(sg_geom_off + 11 * ccomps * dcomps);

            auto g_y_0_yyzz = cbuffer.data(sg_geom_off + 12 * ccomps * dcomps);

            auto g_y_0_yzzz = cbuffer.data(sg_geom_off + 13 * ccomps * dcomps);

            auto g_y_0_zzzz = cbuffer.data(sg_geom_off + 14 * ccomps * dcomps);
            
            sg_geom_off += 15 * ccomps * dcomps;
            
            auto g_z_0_xxxx = cbuffer.data(sg_geom_off + 0 * ccomps * dcomps);

            auto g_z_0_xxxy = cbuffer.data(sg_geom_off + 1 * ccomps * dcomps);

            auto g_z_0_xxxz = cbuffer.data(sg_geom_off + 2 * ccomps * dcomps);

            auto g_z_0_xxyy = cbuffer.data(sg_geom_off + 3 * ccomps * dcomps);

            auto g_z_0_xxyz = cbuffer.data(sg_geom_off + 4 * ccomps * dcomps);

            auto g_z_0_xxzz = cbuffer.data(sg_geom_off + 5 * ccomps * dcomps);

            auto g_z_0_xyyy = cbuffer.data(sg_geom_off + 6 * ccomps * dcomps);

            auto g_z_0_xyyz = cbuffer.data(sg_geom_off + 7 * ccomps * dcomps);

            auto g_z_0_xyzz = cbuffer.data(sg_geom_off + 8 * ccomps * dcomps);

            auto g_z_0_xzzz = cbuffer.data(sg_geom_off + 9 * ccomps * dcomps);

            auto g_z_0_yyyy = cbuffer.data(sg_geom_off + 10 * ccomps * dcomps);

            auto g_z_0_yyyz = cbuffer.data(sg_geom_off + 11 * ccomps * dcomps);

            auto g_z_0_yyzz = cbuffer.data(sg_geom_off + 12 * ccomps * dcomps);

            auto g_z_0_yzzz = cbuffer.data(sg_geom_off + 13 * ccomps * dcomps);

            auto g_z_0_zzzz = cbuffer.data(sg_geom_off + 14 * ccomps * dcomps);
            
            /// Set up components of auxilary buffer : SFSS
            
            auto sf_off = idx_sfxx + i * dcomps + j;

            auto g_0_xxx = cbuffer.data(sf_off + 0 * ccomps * dcomps);

            auto g_0_xxy = cbuffer.data(sf_off + 1 * ccomps * dcomps);

            auto g_0_xxz = cbuffer.data(sf_off + 2 * ccomps * dcomps);

            auto g_0_xyy = cbuffer.data(sf_off + 3 * ccomps * dcomps);

            auto g_0_xyz = cbuffer.data(sf_off + 4 * ccomps * dcomps);

            auto g_0_xzz = cbuffer.data(sf_off + 5 * ccomps * dcomps);

            auto g_0_yyy = cbuffer.data(sf_off + 6 * ccomps * dcomps);

            auto g_0_yyz = cbuffer.data(sf_off + 7 * ccomps * dcomps);

            auto g_0_yzz = cbuffer.data(sf_off + 8 * ccomps * dcomps);

            auto g_0_zzz = cbuffer.data(sf_off + 9 * ccomps * dcomps);
            
            /// set up bra offset for contr_buffer_pfxx

            auto pf_geom_off = idx_geom_pfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxx = cbuffer.data(pf_geom_off + 0 * ccomps * dcomps);

            auto g_x_x_xxy = cbuffer.data(pf_geom_off + 1 * ccomps * dcomps);

            auto g_x_x_xxz = cbuffer.data(pf_geom_off + 2 * ccomps * dcomps);

            auto g_x_x_xyy = cbuffer.data(pf_geom_off + 3 * ccomps * dcomps);

            auto g_x_x_xyz = cbuffer.data(pf_geom_off + 4 * ccomps * dcomps);

            auto g_x_x_xzz = cbuffer.data(pf_geom_off + 5 * ccomps * dcomps);

            auto g_x_x_yyy = cbuffer.data(pf_geom_off + 6 * ccomps * dcomps);

            auto g_x_x_yyz = cbuffer.data(pf_geom_off + 7 * ccomps * dcomps);

            auto g_x_x_yzz = cbuffer.data(pf_geom_off + 8 * ccomps * dcomps);

            auto g_x_x_zzz = cbuffer.data(pf_geom_off + 9 * ccomps * dcomps);

#pragma omp simd aligned(g_x_0_xxx,      \
                             g_x_0_xxxx, \
                             g_x_0_xxxy, \
                             g_x_0_xxxz, \
                             g_x_0_xxy,  \
                             g_x_0_xxyy, \
                             g_x_0_xxyz, \
                             g_x_0_xxz,  \
                             g_x_0_xxzz, \
                             g_x_0_xyy,  \
                             g_x_0_xyyy, \
                             g_x_0_xyyz, \
                             g_x_0_xyz,  \
                             g_x_0_xyzz, \
                             g_x_0_xzz,  \
                             g_x_0_xzzz, \
                             g_x_0_yyy,  \
                             g_x_0_yyz,  \
                             g_x_0_yzz,  \
                             g_x_0_zzz,  \
                             g_0_xxx,    \
                             g_0_xxy,    \
                             g_0_xxz,    \
                             g_0_xyy,    \
                             g_0_xyz,    \
                             g_0_xzz,    \
                             g_0_yyy,    \
                             g_0_yyz,    \
                             g_0_yzz,    \
                             g_0_zzz,    \
                             g_x_x_xxx,  \
                             g_x_x_xxy,  \
                             g_x_x_xxz,  \
                             g_x_x_xyy,  \
                             g_x_x_xyz,  \
                             g_x_x_xzz,  \
                             g_x_x_yyy,  \
                             g_x_x_yyz,  \
                             g_x_x_yzz,  \
                             g_x_x_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxx[k] = -g_x_0_xxx[k] * ab_x + g_x_0_xxxx[k] - g_0_xxx[k];

                g_x_x_xxy[k] = -g_x_0_xxy[k] * ab_x + g_x_0_xxxy[k] - g_0_xxy[k];

                g_x_x_xxz[k] = -g_x_0_xxz[k] * ab_x + g_x_0_xxxz[k] - g_0_xxz[k];

                g_x_x_xyy[k] = -g_x_0_xyy[k] * ab_x + g_x_0_xxyy[k] - g_0_xyy[k];

                g_x_x_xyz[k] = -g_x_0_xyz[k] * ab_x + g_x_0_xxyz[k] - g_0_xyz[k];

                g_x_x_xzz[k] = -g_x_0_xzz[k] * ab_x + g_x_0_xxzz[k] - g_0_xzz[k];

                g_x_x_yyy[k] = -g_x_0_yyy[k] * ab_x + g_x_0_xyyy[k] - g_0_yyy[k];

                g_x_x_yyz[k] = -g_x_0_yyz[k] * ab_x + g_x_0_xyyz[k] - g_0_yyz[k];

                g_x_x_yzz[k] = -g_x_0_yzz[k] * ab_x + g_x_0_xyzz[k] - g_0_yzz[k];

                g_x_x_zzz[k] = -g_x_0_zzz[k] * ab_x + g_x_0_xzzz[k] - g_0_zzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxx = cbuffer.data(pf_geom_off + 10 * ccomps * dcomps);

            auto g_x_y_xxy = cbuffer.data(pf_geom_off + 11 * ccomps * dcomps);

            auto g_x_y_xxz = cbuffer.data(pf_geom_off + 12 * ccomps * dcomps);

            auto g_x_y_xyy = cbuffer.data(pf_geom_off + 13 * ccomps * dcomps);

            auto g_x_y_xyz = cbuffer.data(pf_geom_off + 14 * ccomps * dcomps);

            auto g_x_y_xzz = cbuffer.data(pf_geom_off + 15 * ccomps * dcomps);

            auto g_x_y_yyy = cbuffer.data(pf_geom_off + 16 * ccomps * dcomps);

            auto g_x_y_yyz = cbuffer.data(pf_geom_off + 17 * ccomps * dcomps);

            auto g_x_y_yzz = cbuffer.data(pf_geom_off + 18 * ccomps * dcomps);

            auto g_x_y_zzz = cbuffer.data(pf_geom_off + 19 * ccomps * dcomps);

#pragma omp simd aligned(g_x_0_xxx,      \
                             g_x_0_xxxy, \
                             g_x_0_xxy,  \
                             g_x_0_xxyy, \
                             g_x_0_xxyz, \
                             g_x_0_xxz,  \
                             g_x_0_xyy,  \
                             g_x_0_xyyy, \
                             g_x_0_xyyz, \
                             g_x_0_xyz,  \
                             g_x_0_xyzz, \
                             g_x_0_xzz,  \
                             g_x_0_yyy,  \
                             g_x_0_yyyy, \
                             g_x_0_yyyz, \
                             g_x_0_yyz,  \
                             g_x_0_yyzz, \
                             g_x_0_yzz,  \
                             g_x_0_yzzz, \
                             g_x_0_zzz,  \
                             g_x_y_xxx,  \
                             g_x_y_xxy,  \
                             g_x_y_xxz,  \
                             g_x_y_xyy,  \
                             g_x_y_xyz,  \
                             g_x_y_xzz,  \
                             g_x_y_yyy,  \
                             g_x_y_yyz,  \
                             g_x_y_yzz,  \
                             g_x_y_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxx[k] = -g_x_0_xxx[k] * ab_y + g_x_0_xxxy[k];

                g_x_y_xxy[k] = -g_x_0_xxy[k] * ab_y + g_x_0_xxyy[k];

                g_x_y_xxz[k] = -g_x_0_xxz[k] * ab_y + g_x_0_xxyz[k];

                g_x_y_xyy[k] = -g_x_0_xyy[k] * ab_y + g_x_0_xyyy[k];

                g_x_y_xyz[k] = -g_x_0_xyz[k] * ab_y + g_x_0_xyyz[k];

                g_x_y_xzz[k] = -g_x_0_xzz[k] * ab_y + g_x_0_xyzz[k];

                g_x_y_yyy[k] = -g_x_0_yyy[k] * ab_y + g_x_0_yyyy[k];

                g_x_y_yyz[k] = -g_x_0_yyz[k] * ab_y + g_x_0_yyyz[k];

                g_x_y_yzz[k] = -g_x_0_yzz[k] * ab_y + g_x_0_yyzz[k];

                g_x_y_zzz[k] = -g_x_0_zzz[k] * ab_y + g_x_0_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxx = cbuffer.data(pf_geom_off + 20 * ccomps * dcomps);

            auto g_x_z_xxy = cbuffer.data(pf_geom_off + 21 * ccomps * dcomps);

            auto g_x_z_xxz = cbuffer.data(pf_geom_off + 22 * ccomps * dcomps);

            auto g_x_z_xyy = cbuffer.data(pf_geom_off + 23 * ccomps * dcomps);

            auto g_x_z_xyz = cbuffer.data(pf_geom_off + 24 * ccomps * dcomps);

            auto g_x_z_xzz = cbuffer.data(pf_geom_off + 25 * ccomps * dcomps);

            auto g_x_z_yyy = cbuffer.data(pf_geom_off + 26 * ccomps * dcomps);

            auto g_x_z_yyz = cbuffer.data(pf_geom_off + 27 * ccomps * dcomps);

            auto g_x_z_yzz = cbuffer.data(pf_geom_off + 28 * ccomps * dcomps);

            auto g_x_z_zzz = cbuffer.data(pf_geom_off + 29 * ccomps * dcomps);

#pragma omp simd aligned(g_x_0_xxx,      \
                             g_x_0_xxxz, \
                             g_x_0_xxy,  \
                             g_x_0_xxyz, \
                             g_x_0_xxz,  \
                             g_x_0_xxzz, \
                             g_x_0_xyy,  \
                             g_x_0_xyyz, \
                             g_x_0_xyz,  \
                             g_x_0_xyzz, \
                             g_x_0_xzz,  \
                             g_x_0_xzzz, \
                             g_x_0_yyy,  \
                             g_x_0_yyyz, \
                             g_x_0_yyz,  \
                             g_x_0_yyzz, \
                             g_x_0_yzz,  \
                             g_x_0_yzzz, \
                             g_x_0_zzz,  \
                             g_x_0_zzzz, \
                             g_x_z_xxx,  \
                             g_x_z_xxy,  \
                             g_x_z_xxz,  \
                             g_x_z_xyy,  \
                             g_x_z_xyz,  \
                             g_x_z_xzz,  \
                             g_x_z_yyy,  \
                             g_x_z_yyz,  \
                             g_x_z_yzz,  \
                             g_x_z_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxx[k] = -g_x_0_xxx[k] * ab_z + g_x_0_xxxz[k];

                g_x_z_xxy[k] = -g_x_0_xxy[k] * ab_z + g_x_0_xxyz[k];

                g_x_z_xxz[k] = -g_x_0_xxz[k] * ab_z + g_x_0_xxzz[k];

                g_x_z_xyy[k] = -g_x_0_xyy[k] * ab_z + g_x_0_xyyz[k];

                g_x_z_xyz[k] = -g_x_0_xyz[k] * ab_z + g_x_0_xyzz[k];

                g_x_z_xzz[k] = -g_x_0_xzz[k] * ab_z + g_x_0_xzzz[k];

                g_x_z_yyy[k] = -g_x_0_yyy[k] * ab_z + g_x_0_yyyz[k];

                g_x_z_yyz[k] = -g_x_0_yyz[k] * ab_z + g_x_0_yyzz[k];

                g_x_z_yzz[k] = -g_x_0_yzz[k] * ab_z + g_x_0_yzzz[k];

                g_x_z_zzz[k] = -g_x_0_zzz[k] * ab_z + g_x_0_zzzz[k];
            }
            
            // updated offset
            
            pf_geom_off += 30 * ccomps * dcomps;
            
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxx = cbuffer.data(pf_geom_off + 0 * ccomps * dcomps);

            auto g_y_x_xxy = cbuffer.data(pf_geom_off + 1 * ccomps * dcomps);

            auto g_y_x_xxz = cbuffer.data(pf_geom_off + 2 * ccomps * dcomps);

            auto g_y_x_xyy = cbuffer.data(pf_geom_off + 3 * ccomps * dcomps);

            auto g_y_x_xyz = cbuffer.data(pf_geom_off + 4 * ccomps * dcomps);

            auto g_y_x_xzz = cbuffer.data(pf_geom_off + 5 * ccomps * dcomps);

            auto g_y_x_yyy = cbuffer.data(pf_geom_off + 6 * ccomps * dcomps);

            auto g_y_x_yyz = cbuffer.data(pf_geom_off + 7 * ccomps * dcomps);

            auto g_y_x_yzz = cbuffer.data(pf_geom_off + 8 * ccomps * dcomps);

            auto g_y_x_zzz = cbuffer.data(pf_geom_off + 9 * ccomps * dcomps);

#pragma omp simd aligned(g_y_0_xxx,      \
                             g_y_0_xxxx, \
                             g_y_0_xxxy, \
                             g_y_0_xxxz, \
                             g_y_0_xxy,  \
                             g_y_0_xxyy, \
                             g_y_0_xxyz, \
                             g_y_0_xxz,  \
                             g_y_0_xxzz, \
                             g_y_0_xyy,  \
                             g_y_0_xyyy, \
                             g_y_0_xyyz, \
                             g_y_0_xyz,  \
                             g_y_0_xyzz, \
                             g_y_0_xzz,  \
                             g_y_0_xzzz, \
                             g_y_0_yyy,  \
                             g_y_0_yyz,  \
                             g_y_0_yzz,  \
                             g_y_0_zzz,  \
                             g_y_x_xxx,  \
                             g_y_x_xxy,  \
                             g_y_x_xxz,  \
                             g_y_x_xyy,  \
                             g_y_x_xyz,  \
                             g_y_x_xzz,  \
                             g_y_x_yyy,  \
                             g_y_x_yyz,  \
                             g_y_x_yzz,  \
                             g_y_x_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxx[k] = -g_y_0_xxx[k] * ab_x + g_y_0_xxxx[k];

                g_y_x_xxy[k] = -g_y_0_xxy[k] * ab_x + g_y_0_xxxy[k];

                g_y_x_xxz[k] = -g_y_0_xxz[k] * ab_x + g_y_0_xxxz[k];

                g_y_x_xyy[k] = -g_y_0_xyy[k] * ab_x + g_y_0_xxyy[k];

                g_y_x_xyz[k] = -g_y_0_xyz[k] * ab_x + g_y_0_xxyz[k];

                g_y_x_xzz[k] = -g_y_0_xzz[k] * ab_x + g_y_0_xxzz[k];

                g_y_x_yyy[k] = -g_y_0_yyy[k] * ab_x + g_y_0_xyyy[k];

                g_y_x_yyz[k] = -g_y_0_yyz[k] * ab_x + g_y_0_xyyz[k];

                g_y_x_yzz[k] = -g_y_0_yzz[k] * ab_x + g_y_0_xyzz[k];

                g_y_x_zzz[k] = -g_y_0_zzz[k] * ab_x + g_y_0_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxx = cbuffer.data(pf_geom_off + 10 * ccomps * dcomps);

            auto g_y_y_xxy = cbuffer.data(pf_geom_off + 11 * ccomps * dcomps);

            auto g_y_y_xxz = cbuffer.data(pf_geom_off + 12 * ccomps * dcomps);

            auto g_y_y_xyy = cbuffer.data(pf_geom_off + 13 * ccomps * dcomps);

            auto g_y_y_xyz = cbuffer.data(pf_geom_off + 14 * ccomps * dcomps);

            auto g_y_y_xzz = cbuffer.data(pf_geom_off + 15 * ccomps * dcomps);

            auto g_y_y_yyy = cbuffer.data(pf_geom_off + 16 * ccomps * dcomps);

            auto g_y_y_yyz = cbuffer.data(pf_geom_off + 17 * ccomps * dcomps);

            auto g_y_y_yzz = cbuffer.data(pf_geom_off + 18 * ccomps * dcomps);

            auto g_y_y_zzz = cbuffer.data(pf_geom_off + 19 * ccomps * dcomps);

#pragma omp simd aligned(g_y_0_xxx,      \
                             g_y_0_xxxy, \
                             g_y_0_xxy,  \
                             g_y_0_xxyy, \
                             g_y_0_xxyz, \
                             g_y_0_xxz,  \
                             g_y_0_xyy,  \
                             g_y_0_xyyy, \
                             g_y_0_xyyz, \
                             g_y_0_xyz,  \
                             g_y_0_xyzz, \
                             g_y_0_xzz,  \
                             g_y_0_yyy,  \
                             g_y_0_yyyy, \
                             g_y_0_yyyz, \
                             g_y_0_yyz,  \
                             g_y_0_yyzz, \
                             g_y_0_yzz,  \
                             g_y_0_yzzz, \
                             g_y_0_zzz,  \
                             g_0_xxx,    \
                             g_0_xxy,    \
                             g_0_xxz,    \
                             g_0_xyy,    \
                             g_0_xyz,    \
                             g_0_xzz,    \
                             g_0_yyy,    \
                             g_0_yyz,    \
                             g_0_yzz,    \
                             g_0_zzz,    \
                             g_y_y_xxx,  \
                             g_y_y_xxy,  \
                             g_y_y_xxz,  \
                             g_y_y_xyy,  \
                             g_y_y_xyz,  \
                             g_y_y_xzz,  \
                             g_y_y_yyy,  \
                             g_y_y_yyz,  \
                             g_y_y_yzz,  \
                             g_y_y_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxx[k] = -g_y_0_xxx[k] * ab_y + g_y_0_xxxy[k] - g_0_xxx[k];

                g_y_y_xxy[k] = -g_y_0_xxy[k] * ab_y + g_y_0_xxyy[k] - g_0_xxy[k];

                g_y_y_xxz[k] = -g_y_0_xxz[k] * ab_y + g_y_0_xxyz[k] - g_0_xxz[k];

                g_y_y_xyy[k] = -g_y_0_xyy[k] * ab_y + g_y_0_xyyy[k] - g_0_xyy[k];

                g_y_y_xyz[k] = -g_y_0_xyz[k] * ab_y + g_y_0_xyyz[k] - g_0_xyz[k];

                g_y_y_xzz[k] = -g_y_0_xzz[k] * ab_y + g_y_0_xyzz[k] - g_0_xzz[k];

                g_y_y_yyy[k] = -g_y_0_yyy[k] * ab_y + g_y_0_yyyy[k] - g_0_yyy[k];

                g_y_y_yyz[k] = -g_y_0_yyz[k] * ab_y + g_y_0_yyyz[k] - g_0_yyz[k];

                g_y_y_yzz[k] = -g_y_0_yzz[k] * ab_y + g_y_0_yyzz[k] - g_0_yzz[k];

                g_y_y_zzz[k] = -g_y_0_zzz[k] * ab_y + g_y_0_yzzz[k] - g_0_zzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxx = cbuffer.data(pf_geom_off + 20 * ccomps * dcomps);

            auto g_y_z_xxy = cbuffer.data(pf_geom_off + 21 * ccomps * dcomps);

            auto g_y_z_xxz = cbuffer.data(pf_geom_off + 22 * ccomps * dcomps);

            auto g_y_z_xyy = cbuffer.data(pf_geom_off + 23 * ccomps * dcomps);

            auto g_y_z_xyz = cbuffer.data(pf_geom_off + 24 * ccomps * dcomps);

            auto g_y_z_xzz = cbuffer.data(pf_geom_off + 25 * ccomps * dcomps);

            auto g_y_z_yyy = cbuffer.data(pf_geom_off + 26 * ccomps * dcomps);

            auto g_y_z_yyz = cbuffer.data(pf_geom_off + 27 * ccomps * dcomps);

            auto g_y_z_yzz = cbuffer.data(pf_geom_off + 28 * ccomps * dcomps);

            auto g_y_z_zzz = cbuffer.data(pf_geom_off + 29 * ccomps * dcomps);

#pragma omp simd aligned(g_y_0_xxx,      \
                             g_y_0_xxxz, \
                             g_y_0_xxy,  \
                             g_y_0_xxyz, \
                             g_y_0_xxz,  \
                             g_y_0_xxzz, \
                             g_y_0_xyy,  \
                             g_y_0_xyyz, \
                             g_y_0_xyz,  \
                             g_y_0_xyzz, \
                             g_y_0_xzz,  \
                             g_y_0_xzzz, \
                             g_y_0_yyy,  \
                             g_y_0_yyyz, \
                             g_y_0_yyz,  \
                             g_y_0_yyzz, \
                             g_y_0_yzz,  \
                             g_y_0_yzzz, \
                             g_y_0_zzz,  \
                             g_y_0_zzzz, \
                             g_y_z_xxx,  \
                             g_y_z_xxy,  \
                             g_y_z_xxz,  \
                             g_y_z_xyy,  \
                             g_y_z_xyz,  \
                             g_y_z_xzz,  \
                             g_y_z_yyy,  \
                             g_y_z_yyz,  \
                             g_y_z_yzz,  \
                             g_y_z_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxx[k] = -g_y_0_xxx[k] * ab_z + g_y_0_xxxz[k];

                g_y_z_xxy[k] = -g_y_0_xxy[k] * ab_z + g_y_0_xxyz[k];

                g_y_z_xxz[k] = -g_y_0_xxz[k] * ab_z + g_y_0_xxzz[k];

                g_y_z_xyy[k] = -g_y_0_xyy[k] * ab_z + g_y_0_xyyz[k];

                g_y_z_xyz[k] = -g_y_0_xyz[k] * ab_z + g_y_0_xyzz[k];

                g_y_z_xzz[k] = -g_y_0_xzz[k] * ab_z + g_y_0_xzzz[k];

                g_y_z_yyy[k] = -g_y_0_yyy[k] * ab_z + g_y_0_yyyz[k];

                g_y_z_yyz[k] = -g_y_0_yyz[k] * ab_z + g_y_0_yyzz[k];

                g_y_z_yzz[k] = -g_y_0_yzz[k] * ab_z + g_y_0_yzzz[k];

                g_y_z_zzz[k] = -g_y_0_zzz[k] * ab_z + g_y_0_zzzz[k];
            }
            
            // updated offset
            
            pf_geom_off += 30 * ccomps * dcomps;
            
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxx = cbuffer.data(pf_geom_off + 0 * ccomps * dcomps);

            auto g_z_x_xxy = cbuffer.data(pf_geom_off + 1 * ccomps * dcomps);

            auto g_z_x_xxz = cbuffer.data(pf_geom_off + 2 * ccomps * dcomps);

            auto g_z_x_xyy = cbuffer.data(pf_geom_off + 3 * ccomps * dcomps);

            auto g_z_x_xyz = cbuffer.data(pf_geom_off + 4 * ccomps * dcomps);

            auto g_z_x_xzz = cbuffer.data(pf_geom_off + 5 * ccomps * dcomps);

            auto g_z_x_yyy = cbuffer.data(pf_geom_off + 6 * ccomps * dcomps);

            auto g_z_x_yyz = cbuffer.data(pf_geom_off + 7 * ccomps * dcomps);

            auto g_z_x_yzz = cbuffer.data(pf_geom_off + 8 * ccomps * dcomps);

            auto g_z_x_zzz = cbuffer.data(pf_geom_off + 9 * ccomps * dcomps);

#pragma omp simd aligned(g_z_0_xxx,      \
                             g_z_0_xxxx, \
                             g_z_0_xxxy, \
                             g_z_0_xxxz, \
                             g_z_0_xxy,  \
                             g_z_0_xxyy, \
                             g_z_0_xxyz, \
                             g_z_0_xxz,  \
                             g_z_0_xxzz, \
                             g_z_0_xyy,  \
                             g_z_0_xyyy, \
                             g_z_0_xyyz, \
                             g_z_0_xyz,  \
                             g_z_0_xyzz, \
                             g_z_0_xzz,  \
                             g_z_0_xzzz, \
                             g_z_0_yyy,  \
                             g_z_0_yyz,  \
                             g_z_0_yzz,  \
                             g_z_0_zzz,  \
                             g_z_x_xxx,  \
                             g_z_x_xxy,  \
                             g_z_x_xxz,  \
                             g_z_x_xyy,  \
                             g_z_x_xyz,  \
                             g_z_x_xzz,  \
                             g_z_x_yyy,  \
                             g_z_x_yyz,  \
                             g_z_x_yzz,  \
                             g_z_x_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxx[k] = -g_z_0_xxx[k] * ab_x + g_z_0_xxxx[k];

                g_z_x_xxy[k] = -g_z_0_xxy[k] * ab_x + g_z_0_xxxy[k];

                g_z_x_xxz[k] = -g_z_0_xxz[k] * ab_x + g_z_0_xxxz[k];

                g_z_x_xyy[k] = -g_z_0_xyy[k] * ab_x + g_z_0_xxyy[k];

                g_z_x_xyz[k] = -g_z_0_xyz[k] * ab_x + g_z_0_xxyz[k];

                g_z_x_xzz[k] = -g_z_0_xzz[k] * ab_x + g_z_0_xxzz[k];

                g_z_x_yyy[k] = -g_z_0_yyy[k] * ab_x + g_z_0_xyyy[k];

                g_z_x_yyz[k] = -g_z_0_yyz[k] * ab_x + g_z_0_xyyz[k];

                g_z_x_yzz[k] = -g_z_0_yzz[k] * ab_x + g_z_0_xyzz[k];

                g_z_x_zzz[k] = -g_z_0_zzz[k] * ab_x + g_z_0_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxx = cbuffer.data(pf_geom_off + 10 * ccomps * dcomps);

            auto g_z_y_xxy = cbuffer.data(pf_geom_off + 11 * ccomps * dcomps);

            auto g_z_y_xxz = cbuffer.data(pf_geom_off + 12 * ccomps * dcomps);

            auto g_z_y_xyy = cbuffer.data(pf_geom_off + 13 * ccomps * dcomps);

            auto g_z_y_xyz = cbuffer.data(pf_geom_off + 14 * ccomps * dcomps);

            auto g_z_y_xzz = cbuffer.data(pf_geom_off + 15 * ccomps * dcomps);

            auto g_z_y_yyy = cbuffer.data(pf_geom_off + 16 * ccomps * dcomps);

            auto g_z_y_yyz = cbuffer.data(pf_geom_off + 17 * ccomps * dcomps);

            auto g_z_y_yzz = cbuffer.data(pf_geom_off + 18 * ccomps * dcomps);

            auto g_z_y_zzz = cbuffer.data(pf_geom_off + 19 * ccomps * dcomps);

#pragma omp simd aligned(g_z_0_xxx,      \
                             g_z_0_xxxy, \
                             g_z_0_xxy,  \
                             g_z_0_xxyy, \
                             g_z_0_xxyz, \
                             g_z_0_xxz,  \
                             g_z_0_xyy,  \
                             g_z_0_xyyy, \
                             g_z_0_xyyz, \
                             g_z_0_xyz,  \
                             g_z_0_xyzz, \
                             g_z_0_xzz,  \
                             g_z_0_yyy,  \
                             g_z_0_yyyy, \
                             g_z_0_yyyz, \
                             g_z_0_yyz,  \
                             g_z_0_yyzz, \
                             g_z_0_yzz,  \
                             g_z_0_yzzz, \
                             g_z_0_zzz,  \
                             g_z_y_xxx,  \
                             g_z_y_xxy,  \
                             g_z_y_xxz,  \
                             g_z_y_xyy,  \
                             g_z_y_xyz,  \
                             g_z_y_xzz,  \
                             g_z_y_yyy,  \
                             g_z_y_yyz,  \
                             g_z_y_yzz,  \
                             g_z_y_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxx[k] = -g_z_0_xxx[k] * ab_y + g_z_0_xxxy[k];

                g_z_y_xxy[k] = -g_z_0_xxy[k] * ab_y + g_z_0_xxyy[k];

                g_z_y_xxz[k] = -g_z_0_xxz[k] * ab_y + g_z_0_xxyz[k];

                g_z_y_xyy[k] = -g_z_0_xyy[k] * ab_y + g_z_0_xyyy[k];

                g_z_y_xyz[k] = -g_z_0_xyz[k] * ab_y + g_z_0_xyyz[k];

                g_z_y_xzz[k] = -g_z_0_xzz[k] * ab_y + g_z_0_xyzz[k];

                g_z_y_yyy[k] = -g_z_0_yyy[k] * ab_y + g_z_0_yyyy[k];

                g_z_y_yyz[k] = -g_z_0_yyz[k] * ab_y + g_z_0_yyyz[k];

                g_z_y_yzz[k] = -g_z_0_yzz[k] * ab_y + g_z_0_yyzz[k];

                g_z_y_zzz[k] = -g_z_0_zzz[k] * ab_y + g_z_0_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxx = cbuffer.data(pf_geom_off + 20 * ccomps * dcomps);

            auto g_z_z_xxy = cbuffer.data(pf_geom_off + 21 * ccomps * dcomps);

            auto g_z_z_xxz = cbuffer.data(pf_geom_off + 22 * ccomps * dcomps);

            auto g_z_z_xyy = cbuffer.data(pf_geom_off + 23 * ccomps * dcomps);

            auto g_z_z_xyz = cbuffer.data(pf_geom_off + 24 * ccomps * dcomps);

            auto g_z_z_xzz = cbuffer.data(pf_geom_off + 25 * ccomps * dcomps);

            auto g_z_z_yyy = cbuffer.data(pf_geom_off + 26 * ccomps * dcomps);

            auto g_z_z_yyz = cbuffer.data(pf_geom_off + 27 * ccomps * dcomps);

            auto g_z_z_yzz = cbuffer.data(pf_geom_off + 28 * ccomps * dcomps);

            auto g_z_z_zzz = cbuffer.data(pf_geom_off + 29 * ccomps * dcomps);

#pragma omp simd aligned(g_z_0_xxx,      \
                             g_z_0_xxxz, \
                             g_z_0_xxy,  \
                             g_z_0_xxyz, \
                             g_z_0_xxz,  \
                             g_z_0_xxzz, \
                             g_z_0_xyy,  \
                             g_z_0_xyyz, \
                             g_z_0_xyz,  \
                             g_z_0_xyzz, \
                             g_z_0_xzz,  \
                             g_z_0_xzzz, \
                             g_z_0_yyy,  \
                             g_z_0_yyyz, \
                             g_z_0_yyz,  \
                             g_z_0_yyzz, \
                             g_z_0_yzz,  \
                             g_z_0_yzzz, \
                             g_z_0_zzz,  \
                             g_z_0_zzzz, \
                             g_0_xxx,    \
                             g_0_xxy,    \
                             g_0_xxz,    \
                             g_0_xyy,    \
                             g_0_xyz,    \
                             g_0_xzz,    \
                             g_0_yyy,    \
                             g_0_yyz,    \
                             g_0_yzz,    \
                             g_0_zzz,    \
                             g_z_z_xxx,  \
                             g_z_z_xxy,  \
                             g_z_z_xxz,  \
                             g_z_z_xyy,  \
                             g_z_z_xyz,  \
                             g_z_z_xzz,  \
                             g_z_z_yyy,  \
                             g_z_z_yyz,  \
                             g_z_z_yzz,  \
                             g_z_z_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxx[k] = -g_z_0_xxx[k] * ab_z + g_z_0_xxxz[k] - g_0_xxx[k];

                g_z_z_xxy[k] = -g_z_0_xxy[k] * ab_z + g_z_0_xxyz[k] - g_0_xxy[k];

                g_z_z_xxz[k] = -g_z_0_xxz[k] * ab_z + g_z_0_xxzz[k] - g_0_xxz[k];

                g_z_z_xyy[k] = -g_z_0_xyy[k] * ab_z + g_z_0_xyyz[k] - g_0_xyy[k];

                g_z_z_xyz[k] = -g_z_0_xyz[k] * ab_z + g_z_0_xyzz[k] - g_0_xyz[k];

                g_z_z_xzz[k] = -g_z_0_xzz[k] * ab_z + g_z_0_xzzz[k] - g_0_xzz[k];

                g_z_z_yyy[k] = -g_z_0_yyy[k] * ab_z + g_z_0_yyyz[k] - g_0_yyy[k];

                g_z_z_yyz[k] = -g_z_0_yyz[k] * ab_z + g_z_0_yyzz[k] - g_0_yyz[k];

                g_z_z_yzz[k] = -g_z_0_yzz[k] * ab_z + g_z_0_yzzz[k] - g_0_yzz[k];

                g_z_z_zzz[k] = -g_z_0_zzz[k] * ab_z + g_z_0_zzzz[k] - g_0_zzz[k];
            }
        }
    }
}

}  // namespace erirec
