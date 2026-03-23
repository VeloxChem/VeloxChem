#include "LocalCorePotentialPrimRecPG.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_pg(CSimdArray<double>& pbuffer, 
                                  const size_t idx_pg,
                                  const size_t idx_sf,
                                  const size_t idx_sg,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx = pbuffer.data(idx_sf);

    auto tg_0_xxy = pbuffer.data(idx_sf + 1);

    auto tg_0_xxz = pbuffer.data(idx_sf + 2);

    auto tg_0_xyy = pbuffer.data(idx_sf + 3);

    auto tg_0_xyz = pbuffer.data(idx_sf + 4);

    auto tg_0_xzz = pbuffer.data(idx_sf + 5);

    auto tg_0_yyy = pbuffer.data(idx_sf + 6);

    auto tg_0_yyz = pbuffer.data(idx_sf + 7);

    auto tg_0_yzz = pbuffer.data(idx_sf + 8);

    auto tg_0_zzz = pbuffer.data(idx_sf + 9);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx = pbuffer.data(idx_sg);

    auto tg_0_xxxy = pbuffer.data(idx_sg + 1);

    auto tg_0_xxxz = pbuffer.data(idx_sg + 2);

    auto tg_0_xxyy = pbuffer.data(idx_sg + 3);

    auto tg_0_xxyz = pbuffer.data(idx_sg + 4);

    auto tg_0_xxzz = pbuffer.data(idx_sg + 5);

    auto tg_0_xyyy = pbuffer.data(idx_sg + 6);

    auto tg_0_xyyz = pbuffer.data(idx_sg + 7);

    auto tg_0_xyzz = pbuffer.data(idx_sg + 8);

    auto tg_0_xzzz = pbuffer.data(idx_sg + 9);

    auto tg_0_yyyy = pbuffer.data(idx_sg + 10);

    auto tg_0_yyyz = pbuffer.data(idx_sg + 11);

    auto tg_0_yyzz = pbuffer.data(idx_sg + 12);

    auto tg_0_yzzz = pbuffer.data(idx_sg + 13);

    auto tg_0_zzzz = pbuffer.data(idx_sg + 14);

    // Set up components of targeted buffer : PG

    auto tg_x_xxxx = pbuffer.data(idx_pg);

    auto tg_x_xxxy = pbuffer.data(idx_pg + 1);

    auto tg_x_xxxz = pbuffer.data(idx_pg + 2);

    auto tg_x_xxyy = pbuffer.data(idx_pg + 3);

    auto tg_x_xxyz = pbuffer.data(idx_pg + 4);

    auto tg_x_xxzz = pbuffer.data(idx_pg + 5);

    auto tg_x_xyyy = pbuffer.data(idx_pg + 6);

    auto tg_x_xyyz = pbuffer.data(idx_pg + 7);

    auto tg_x_xyzz = pbuffer.data(idx_pg + 8);

    auto tg_x_xzzz = pbuffer.data(idx_pg + 9);

    auto tg_x_yyyy = pbuffer.data(idx_pg + 10);

    auto tg_x_yyyz = pbuffer.data(idx_pg + 11);

    auto tg_x_yyzz = pbuffer.data(idx_pg + 12);

    auto tg_x_yzzz = pbuffer.data(idx_pg + 13);

    auto tg_x_zzzz = pbuffer.data(idx_pg + 14);

    auto tg_y_xxxx = pbuffer.data(idx_pg + 15);

    auto tg_y_xxxy = pbuffer.data(idx_pg + 16);

    auto tg_y_xxxz = pbuffer.data(idx_pg + 17);

    auto tg_y_xxyy = pbuffer.data(idx_pg + 18);

    auto tg_y_xxyz = pbuffer.data(idx_pg + 19);

    auto tg_y_xxzz = pbuffer.data(idx_pg + 20);

    auto tg_y_xyyy = pbuffer.data(idx_pg + 21);

    auto tg_y_xyyz = pbuffer.data(idx_pg + 22);

    auto tg_y_xyzz = pbuffer.data(idx_pg + 23);

    auto tg_y_xzzz = pbuffer.data(idx_pg + 24);

    auto tg_y_yyyy = pbuffer.data(idx_pg + 25);

    auto tg_y_yyyz = pbuffer.data(idx_pg + 26);

    auto tg_y_yyzz = pbuffer.data(idx_pg + 27);

    auto tg_y_yzzz = pbuffer.data(idx_pg + 28);

    auto tg_y_zzzz = pbuffer.data(idx_pg + 29);

    auto tg_z_xxxx = pbuffer.data(idx_pg + 30);

    auto tg_z_xxxy = pbuffer.data(idx_pg + 31);

    auto tg_z_xxxz = pbuffer.data(idx_pg + 32);

    auto tg_z_xxyy = pbuffer.data(idx_pg + 33);

    auto tg_z_xxyz = pbuffer.data(idx_pg + 34);

    auto tg_z_xxzz = pbuffer.data(idx_pg + 35);

    auto tg_z_xyyy = pbuffer.data(idx_pg + 36);

    auto tg_z_xyyz = pbuffer.data(idx_pg + 37);

    auto tg_z_xyzz = pbuffer.data(idx_pg + 38);

    auto tg_z_xzzz = pbuffer.data(idx_pg + 39);

    auto tg_z_yyyy = pbuffer.data(idx_pg + 40);

    auto tg_z_yyyz = pbuffer.data(idx_pg + 41);

    auto tg_z_yyzz = pbuffer.data(idx_pg + 42);

    auto tg_z_yzzz = pbuffer.data(idx_pg + 43);

    auto tg_z_zzzz = pbuffer.data(idx_pg + 44);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_xxx, tg_0_xxxx, tg_0_xxxy, tg_0_xxxz, tg_0_xxy, tg_0_xxyy, tg_0_xxyz, tg_0_xxz, tg_0_xxzz, tg_0_xyy, tg_0_xyyy, tg_0_xyyz, tg_0_xyz, tg_0_xyzz, tg_0_xzz, tg_0_xzzz, tg_0_yyy, tg_0_yyyy, tg_0_yyyz, tg_0_yyz, tg_0_yyzz, tg_0_yzz, tg_0_yzzz, tg_0_zzz, tg_0_zzzz, tg_x_xxxx, tg_x_xxxy, tg_x_xxxz, tg_x_xxyy, tg_x_xxyz, tg_x_xxzz, tg_x_xyyy, tg_x_xyyz, tg_x_xyzz, tg_x_xzzz, tg_x_yyyy, tg_x_yyyz, tg_x_yyzz, tg_x_yzzz, tg_x_zzzz, tg_y_xxxx, tg_y_xxxy, tg_y_xxxz, tg_y_xxyy, tg_y_xxyz, tg_y_xxzz, tg_y_xyyy, tg_y_xyyz, tg_y_xyzz, tg_y_xzzz, tg_y_yyyy, tg_y_yyyz, tg_y_yyzz, tg_y_yzzz, tg_y_zzzz, tg_z_xxxx, tg_z_xxxy, tg_z_xxxz, tg_z_xxyy, tg_z_xxyz, tg_z_xxzz, tg_z_xyyy, tg_z_xyyz, tg_z_xyzz, tg_z_xzzz, tg_z_yyyy, tg_z_yyyz, tg_z_yyzz, tg_z_yzzz, tg_z_zzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_x_xxxx[i] = 4.0 * tg_0_xxx[i] * fxi[i] + tg_0_xxxx[i] * ra_x[i];

        tg_x_xxxy[i] = 3.0 * tg_0_xxy[i] * fxi[i] + tg_0_xxxy[i] * ra_x[i];

        tg_x_xxxz[i] = 3.0 * tg_0_xxz[i] * fxi[i] + tg_0_xxxz[i] * ra_x[i];

        tg_x_xxyy[i] = 2.0 * tg_0_xyy[i] * fxi[i] + tg_0_xxyy[i] * ra_x[i];

        tg_x_xxyz[i] = 2.0 * tg_0_xyz[i] * fxi[i] + tg_0_xxyz[i] * ra_x[i];

        tg_x_xxzz[i] = 2.0 * tg_0_xzz[i] * fxi[i] + tg_0_xxzz[i] * ra_x[i];

        tg_x_xyyy[i] = tg_0_yyy[i] * fxi[i] + tg_0_xyyy[i] * ra_x[i];

        tg_x_xyyz[i] = tg_0_yyz[i] * fxi[i] + tg_0_xyyz[i] * ra_x[i];

        tg_x_xyzz[i] = tg_0_yzz[i] * fxi[i] + tg_0_xyzz[i] * ra_x[i];

        tg_x_xzzz[i] = tg_0_zzz[i] * fxi[i] + tg_0_xzzz[i] * ra_x[i];

        tg_x_yyyy[i] = tg_0_yyyy[i] * ra_x[i];

        tg_x_yyyz[i] = tg_0_yyyz[i] * ra_x[i];

        tg_x_yyzz[i] = tg_0_yyzz[i] * ra_x[i];

        tg_x_yzzz[i] = tg_0_yzzz[i] * ra_x[i];

        tg_x_zzzz[i] = tg_0_zzzz[i] * ra_x[i];

        tg_y_xxxx[i] = tg_0_xxxx[i] * ra_y[i];

        tg_y_xxxy[i] = tg_0_xxx[i] * fxi[i] + tg_0_xxxy[i] * ra_y[i];

        tg_y_xxxz[i] = tg_0_xxxz[i] * ra_y[i];

        tg_y_xxyy[i] = 2.0 * tg_0_xxy[i] * fxi[i] + tg_0_xxyy[i] * ra_y[i];

        tg_y_xxyz[i] = tg_0_xxz[i] * fxi[i] + tg_0_xxyz[i] * ra_y[i];

        tg_y_xxzz[i] = tg_0_xxzz[i] * ra_y[i];

        tg_y_xyyy[i] = 3.0 * tg_0_xyy[i] * fxi[i] + tg_0_xyyy[i] * ra_y[i];

        tg_y_xyyz[i] = 2.0 * tg_0_xyz[i] * fxi[i] + tg_0_xyyz[i] * ra_y[i];

        tg_y_xyzz[i] = tg_0_xzz[i] * fxi[i] + tg_0_xyzz[i] * ra_y[i];

        tg_y_xzzz[i] = tg_0_xzzz[i] * ra_y[i];

        tg_y_yyyy[i] = 4.0 * tg_0_yyy[i] * fxi[i] + tg_0_yyyy[i] * ra_y[i];

        tg_y_yyyz[i] = 3.0 * tg_0_yyz[i] * fxi[i] + tg_0_yyyz[i] * ra_y[i];

        tg_y_yyzz[i] = 2.0 * tg_0_yzz[i] * fxi[i] + tg_0_yyzz[i] * ra_y[i];

        tg_y_yzzz[i] = tg_0_zzz[i] * fxi[i] + tg_0_yzzz[i] * ra_y[i];

        tg_y_zzzz[i] = tg_0_zzzz[i] * ra_y[i];

        tg_z_xxxx[i] = tg_0_xxxx[i] * ra_z[i];

        tg_z_xxxy[i] = tg_0_xxxy[i] * ra_z[i];

        tg_z_xxxz[i] = tg_0_xxx[i] * fxi[i] + tg_0_xxxz[i] * ra_z[i];

        tg_z_xxyy[i] = tg_0_xxyy[i] * ra_z[i];

        tg_z_xxyz[i] = tg_0_xxy[i] * fxi[i] + tg_0_xxyz[i] * ra_z[i];

        tg_z_xxzz[i] = 2.0 * tg_0_xxz[i] * fxi[i] + tg_0_xxzz[i] * ra_z[i];

        tg_z_xyyy[i] = tg_0_xyyy[i] * ra_z[i];

        tg_z_xyyz[i] = tg_0_xyy[i] * fxi[i] + tg_0_xyyz[i] * ra_z[i];

        tg_z_xyzz[i] = 2.0 * tg_0_xyz[i] * fxi[i] + tg_0_xyzz[i] * ra_z[i];

        tg_z_xzzz[i] = 3.0 * tg_0_xzz[i] * fxi[i] + tg_0_xzzz[i] * ra_z[i];

        tg_z_yyyy[i] = tg_0_yyyy[i] * ra_z[i];

        tg_z_yyyz[i] = tg_0_yyy[i] * fxi[i] + tg_0_yyyz[i] * ra_z[i];

        tg_z_yyzz[i] = 2.0 * tg_0_yyz[i] * fxi[i] + tg_0_yyzz[i] * ra_z[i];

        tg_z_yzzz[i] = 3.0 * tg_0_yzz[i] * fxi[i] + tg_0_yzzz[i] * ra_z[i];

        tg_z_zzzz[i] = 4.0 * tg_0_zzz[i] * fxi[i] + tg_0_zzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

