#include "LocalCorePotentialPrimRecGS.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_gs(CSimdArray<double>& pbuffer, 
                                  const size_t idx_gs,
                                  const size_t idx_ds,
                                  const size_t idx_fs,
                                  const CSimdArray<double>& factors) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(8);

    auto ra_y = factors.data(9);

    auto ra_z = factors.data(10);

    // Set up inverted 1/2xi

    auto fxi = factors.data(11);

    // Set up components of auxiliary buffer : DS

    auto tg_xx_0 = pbuffer.data(idx_ds);

    auto tg_yy_0 = pbuffer.data(idx_ds + 3);

    auto tg_zz_0 = pbuffer.data(idx_ds + 5);

    // Set up components of auxiliary buffer : FS

    auto tg_xxx_0 = pbuffer.data(idx_fs);

    auto tg_xxz_0 = pbuffer.data(idx_fs + 2);

    auto tg_xyy_0 = pbuffer.data(idx_fs + 3);

    auto tg_xzz_0 = pbuffer.data(idx_fs + 5);

    auto tg_yyy_0 = pbuffer.data(idx_fs + 6);

    auto tg_yyz_0 = pbuffer.data(idx_fs + 7);

    auto tg_yzz_0 = pbuffer.data(idx_fs + 8);

    auto tg_zzz_0 = pbuffer.data(idx_fs + 9);

    // Set up components of targeted buffer : GS

    auto tg_xxxx_0 = pbuffer.data(idx_gs);

    auto tg_xxxy_0 = pbuffer.data(idx_gs + 1);

    auto tg_xxxz_0 = pbuffer.data(idx_gs + 2);

    auto tg_xxyy_0 = pbuffer.data(idx_gs + 3);

    auto tg_xxyz_0 = pbuffer.data(idx_gs + 4);

    auto tg_xxzz_0 = pbuffer.data(idx_gs + 5);

    auto tg_xyyy_0 = pbuffer.data(idx_gs + 6);

    auto tg_xyyz_0 = pbuffer.data(idx_gs + 7);

    auto tg_xyzz_0 = pbuffer.data(idx_gs + 8);

    auto tg_xzzz_0 = pbuffer.data(idx_gs + 9);

    auto tg_yyyy_0 = pbuffer.data(idx_gs + 10);

    auto tg_yyyz_0 = pbuffer.data(idx_gs + 11);

    auto tg_yyzz_0 = pbuffer.data(idx_gs + 12);

    auto tg_yzzz_0 = pbuffer.data(idx_gs + 13);

    auto tg_zzzz_0 = pbuffer.data(idx_gs + 14);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xx_0, tg_xxx_0, tg_xxxx_0, tg_xxxy_0, tg_xxxz_0, tg_xxyy_0, tg_xxyz_0, tg_xxz_0, tg_xxzz_0, tg_xyy_0, tg_xyyy_0, tg_xyyz_0, tg_xyzz_0, tg_xzz_0, tg_xzzz_0, tg_yy_0, tg_yyy_0, tg_yyyy_0, tg_yyyz_0, tg_yyz_0, tg_yyzz_0, tg_yzz_0, tg_yzzz_0, tg_zz_0, tg_zzz_0, tg_zzzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxx_0[i] = 3.0 * tg_xx_0[i] * fxi[i] + tg_xxx_0[i] * ra_x[i];

        tg_xxxy_0[i] = tg_xxx_0[i] * ra_y[i];

        tg_xxxz_0[i] = tg_xxx_0[i] * ra_z[i];

        tg_xxyy_0[i] = tg_yy_0[i] * fxi[i] + tg_xyy_0[i] * ra_x[i];

        tg_xxyz_0[i] = tg_xxz_0[i] * ra_y[i];

        tg_xxzz_0[i] = tg_zz_0[i] * fxi[i] + tg_xzz_0[i] * ra_x[i];

        tg_xyyy_0[i] = tg_yyy_0[i] * ra_x[i];

        tg_xyyz_0[i] = tg_yyz_0[i] * ra_x[i];

        tg_xyzz_0[i] = tg_yzz_0[i] * ra_x[i];

        tg_xzzz_0[i] = tg_zzz_0[i] * ra_x[i];

        tg_yyyy_0[i] = 3.0 * tg_yy_0[i] * fxi[i] + tg_yyy_0[i] * ra_y[i];

        tg_yyyz_0[i] = tg_yyy_0[i] * ra_z[i];

        tg_yyzz_0[i] = tg_zz_0[i] * fxi[i] + tg_yzz_0[i] * ra_y[i];

        tg_yzzz_0[i] = tg_zzz_0[i] * ra_y[i];

        tg_zzzz_0[i] = 3.0 * tg_zz_0[i] * fxi[i] + tg_zzz_0[i] * ra_z[i];
    }
}

} // t2lecp namespace

